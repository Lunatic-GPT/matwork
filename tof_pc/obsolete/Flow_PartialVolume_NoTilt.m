function Flow_PartialVolume_NoTilt(params,save_name)
% modified from Flow_PartialVolume
%% params should contain the following fields:
% roi_vessel: vessel mask
% roi_wm: white matter mask
% voxSize: acquired voxSize; (cm)
% voxSize_interp: voxSize of the images; (cm)

% thk: slice thickness: (cm)
% mag: magnitude images;
% ph: phase images; (value should be in degrees).
% FA: flip angle;
% roiRad: radius for the circular ROI encompassing all pixels to be included in model image fitting (cm).
% bgRad: background ROI out- and inner- radii for calculating tissue signal (cm)
% interp_factor: interpolation factor for PSF convolution.
% VENC: (cm/s)
% T1: T1 for blood and WM tissue, respectively. (s)
% T2: T2 for blood and WM tissue, respectively. (s)
% TR: (s)
% TE: (s)
% negPhase: whether the flow produce a negative phase
% flow_pattern: 'laminar' (default) or 'plug'
% plot_fit: plot fit results
% skip_fit: skip fitting
% v0, rad0: radius and flow velocity in mm and cm/s; initial values for
% fitting; also used for calculating the fitting image when skip_fit is
% true;

    T1=params.T1;  %blood, wm.

    if ~isfield(params,'T2')
     T2=[28.5,23.5]*0.001;
    else
     T2=params.T2; %  %putaman T2* is 17.6 ms
    end
    
    neg_phase=params.negPhase; %false;

bgRad=params.bgRad;

thk=params.thk; % cm

TR = params.TR; % in s

TE = params.TE; %

lambda=params.lambda; % blood to brain water partition coef
%vox_size=[dReadoutFOV,dPhaseFOV]./round([lRO,lPE]*f_nonzero);

voxSize_interp=params.voxSize_interp;


FA=params.FA;
if ~isfield(params,'flow_pattern')
  params.flow_pattern='laminar';
end

[va,sa,s_static]=prep_interp(TR, T1, T2, TE, thk, FA, params.flow_pattern);

%%

roi_vessel=params.roi_vessel;
if max(roi_vessel(:))==1 %not clusterized
 roi_vessel=clusterize2_2d(roi_vessel);
end

roi_wm=params.roi_wm;

nroi=max(roi_vessel(:));
v=zeros(1,nroi);
rad = zeros(1,nroi);

center=zeros(2,nroi);
mbg=zeros(2,nroi);

cd_fit=[];
exitflag=[];
bgroi=[];
data=[];

if neg_phase
    params.ph=-params.ph;
end
              
    do_detrend=true;
for i=1:nroi
    %%
    
    tic;
    roi2=roi_vessel==i;
    if sum(roi2(:))==0
        continue;
    end
    
    bgRad_i=round(bgRad./voxSize_interp);
    
    bg_out=mask_circle(size(roi2),bgRad_i(1),roiCOM(roi2),1);
    bg_in=mask_circle(size(roi2),bgRad_i(2),roiCOM(roi2),1);
         
    bg=bg_out&~bg_in;
    
    
    sl_sel=squeeze(sum(sum(roi2,1),2))>0;
    roi2=roi2(:,:,sl_sel);
    bgroi=bg&roi_wm(:,:,sl_sel)>0&roi_vessel==0;
    
    if sum(bgroi(:))<118 % if outside wm, do not pose restriction.
       bgroi=bg;% 
    end
    
    pos=roiCOM(roi2);
    roiRad=params.roiRad;
    roiRad_i=round(mean(roiRad./voxSize_interp));  % 1 and 2 should be equal
                 
    cropr=pos(1)-roiRad_i:pos(1)+roiRad_i;       
    cropc=pos(2)-roiRad_i:pos(2)+roiRad_i;
            
    cropr_l=pos(1)-bgRad_i(1)-5:pos(1)+bgRad_i(1)+5;       
    cropc_l=pos(2)-bgRad_i(1)-5:pos(2)+bgRad_i(1)+5;
    
    if do_detrend  % do we need to do detrend again?
        if size(params.ph,4)>1 % phase images
         ph=tof_bg_rm(params.ph(:,:,sl_sel,2*i-1:2*i),bgroi,[],[],false);
         mag=tof_bg_rm(params.mag(:,:,sl_sel,2*i-1:2*i),bgroi);
        else % phase contrast images
         ph=tof_bg_rm(params.ph(:,:,sl_sel,:),bgroi,[],[],false);
         mag=tof_bg_rm(params.mag(:,:,sl_sel,:),bgroi);
        end
    else
        mag=params.mag(:,:,sl_sel);
        ph=params.ph(:,:,sl_sel);
    end
    
    
    mbg(:,i)=mean_roi(mag,bgroi>0);               
        s_static=s_static(1);
        
        
        snorm=mean(mbg(:,i))/s_static(1)*lambda;
        
        mag1=mag(:,:,:,1);
        mag2=mag(:,:,:,2);
        
        phbg=0;
        if size(ph,4)==1
           cd=mag2.*exp(1i*(ph-phbg)/180*pi) - mag1;
        else
           cd=mag2.*exp(-1i*ph(:,:,:,2)/180*pi) - mag1.*exp(-1i*ph(:,:,:,1)/180*pi);        
        end
        cd=cd/snorm;
        
        z(1)=mean(cd(roi2));
        
        z(2)=mean(mag1(roi2))/snorm;
    
        ph_crop(:,:,:,i)=ph(cropr_l,cropc_l,:);
        mag1_crop(:,:,:,i)=mag(cropr_l,cropc_l,:,1);
        mag2_crop(:,:,:,i)=mag(cropr_l,cropc_l,:,2);
        cd_crop(:,:,:,i)=cd(cropr_l,cropc_l);
        bgroi_crop(:,:,:,i)=bgroi(cropr_l,cropc_l);
        roi3=mask_circle(size(cd),roiRad_i,pos,1);
        roi_crop(:,:,:,i)=roi3(cropr_l,cropc_l);
          
        cd=cd(cropr,cropc);
        
        
    
        
      
        data(:,:,i)=cd;
            
            
                fit_par=params;
                fit_par.va=va;
                fit_par.sa=sa;
                fit_par.data=data(:,:,i);
                fit_par.roi=roi3(cropr,cropc);
                
                fitfunc2 = @(x) fitfunc(x,fit_par);
                %data(:,:,i),roi3,FOV,vox_size,vox_size./interp_factor,v_calc_pc,sart_calc_pc,VENC,flow_pattern);
                if datenum(version('-date'))>=736580
                    options = optimoptions('lsqnonlin','MaxIter',Inf,'StepTolerance',1e-6,...
                        'MaxFunctionEvaluations',Inf,'FunctionTolerance',0,...
                        'OptimalityTolerance',0,'Display','iter',...
                        'FiniteDifferenceStepSize',0.01);
                else
                    options = optimoptions('lsqnonlin','MaxIter',Inf,'TolX',1e-6,...
                        'MaxFunEvals',Inf,'TolFun',0,'TolPCG',0,'Display','iter','FinDiffRelStep',0.01);
                end
    
    
% x(1:2) in cm
% x(3) in mm 
% x(4) in 10 cm/s

               if isfield(params,'skip_fit') && params.skip_fit
                   res=[0,0,params.rad0,params.v0/10];
               else
                   [res,tmp,tmp2,exitflag(i),output(i)]=lsqnonlin(fitfunc2,[0,0,params.rad0,params.v0/10],[-0.05,-0.05,0.01,0],[0.05,0.05,0.25,0.5],options);
               end
               
               center(:,i) = res(1:2)*10; % in mm
              
            v(i) = res(4)*10;  % in cm/s
            rad(i)=res(3); %% in mm
            
          
            [resid,cd_fit(:,:,i)]=fitfunc2(res);
            
            
            if isfield(params,'plot_fit') && params.plot_fit  %check fit
                disp(res);
                fprintf('diam = %3.2f mm; v = %3.2f cm/s\n',rad(i)*2,v(i)); 
                tmp1=data(:,:,i);
                tmp1(~roi3)=0;
                tmp2=cd_fit(:,:,i);
                tmp2(~roi3)=0;
                
               figure; compare_cd(tmp1,tmp2);
                resid=sos(tmp1(:)-tmp2(:),1).^2;
                disp(resid);

            end
            fprintf('vessel %d/%d: v = %s cm/s; d = %s mm\n',i,nroi,num2str(v(i)),num2str(rad(i)*2));
   fprintf('Time remaining %4.3f s\n',toc*(nroi-i));
end
roi_vessel=params.roi_vessel;
params=rmfield(params,{'mag','ph','roi_vessel'});

if ~exist('output','var')
    output=[];
end
  save(save_name,'v','rad','cd_fit','data','center','params','exitflag','output','bgroi','roi_vessel','ph_crop','mag1_crop','mag2_crop','cd_crop','bgroi_crop','roi_crop');
            
function [v_pc,s_pc,s_static]=prep_interp(TR, T1, T2, TE, thk, FA, flow_pattern)
                
f_prep=['TR_',num2str(TR),'_T1_',num2str(T1),'_T2_',num2str(T2),'_TE_',num2str(TE),'_thk_',num2str(thk),'_FA_',num2str(FA),'_',flow_pattern,'.mat'];
f_prep=strrep(f_prep,' ','_');
for i=1:10
    f_prep=strrep(f_prep,'__','_');
end
    d_prep=fileparts(mfilename('fullpath'));
    
    d_prep=fullfile(d_prep,'prep_interp',f_prep);    
    
if ~exist(d_prep,'file')
    [tmp1,tmp2,v_pc,s_pc,s_static]=Prep_PV_V4TOF2Sur_PC(TR,T1,TE,T2,thk,FA,flow_pattern,'sinc','fl_fq_retroZ_mb');
    save(d_prep,'v_pc','s_pc','s_static');
else
    load(d_prep);
end


 function res=fr(r,vmean,va,sa,rad,venc,flow_pattern)
    if strcmp(flow_pattern,'laminar')
      vi=v_laminarFlow(r/rad,vmean);
    else
        vi=vmean*ones(size(r));
    end
    si=interp1(va,sa,vi);
    res=si.*(exp(1i*vi/venc*pi));
    res(r>rad)=0;
    
    
        function [res,cd_fit]=fitfunc(x,fit_par)
            
% x(1:2) in cm
% x(3) in mm 
% x(4) in 10 cm/s

            %rad,vmean,center,cd,roi,FOV,voxSize,voxSize_interp,va,sa,venc,flow_pattern)
            rad=x(3)*0.1; % to cm
            vmean=x(4)*10; % to cm/s
            center=x(1:2); 
            
            cd = fit_par.data;
            roi= fit_par.roi;
            voxSize=fit_par.voxSize;
            voxSize_interp=fit_par.voxSize_interp;
            venc=fit_par.VENC;
            flow_pattern=fit_par.flow_pattern;
            
            va=fit_par.va;  % cm/s
            sa=fit_par.sa;
            
            FOV=size(cd).*voxSize_interp;

            fr2=@(r) fr(r,vmean,va,sa,rad,venc,flow_pattern);
            im1=circleImage_CartesianKSpace(fr2,center,FOV,voxSize,voxSize_interp,false,fit_par.interp);
            
            fr2=@(r) fr(r,vmean,va,sa,rad,Inf,flow_pattern);
            im0=circleImage_CartesianKSpace(fr2,center,FOV,voxSize,voxSize_interp,false,fit_par.interp);
            
            cd_fit=im1-im0;    
            resid=cd(roi>0)-cd_fit(roi>0);
            res=[real(resid(:));imag(resid(:))];
                
                   
            
            
            
            