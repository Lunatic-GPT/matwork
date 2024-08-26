function Flow_PartialVolume(params,T1,T2,neg_phase)
%meas: raw data name; no suffix
% sub_dir: sub directory data
% rois should be saved in mask_vessels_mid.mat in meas/sub_dir
% white matter mask saved in mask_wm_mid.mat in meas/sub_dir
% phase file should be saved in Phase_%s_mean.mat
% magnitude file should be saved in Mag_%s_mean.mat
% method 1: TOF_PC;
%        2a: Hamilton
%        2b: Hamilton + 2 compartment
%        3a: Cmplex diff
%        3b: Cmplex diff + 2 compartment  (not good; use 2b instead)
%        3c: Cmplex diff: spatial pattern matching
% mask_file: if not given, use the default mask_mid_vessels.mat
% T1 values at 7 T: blood 2.6; wm: 1.2; putaman: 1.7 s.  From Rooney et al.
% T2* values at 7 T: blood 28.5 s; wm 23.5 ms; putaman: 17.6 ms; own
% measurement
% 7/27/2017: changed from Flow_PartialVolume(meas,sub_dir,method,interp_factor,T1,T2,FA_scale,neg_phase,mask_file)
% to Flow_PartialVolume(meas,sub_dir,method,interp_factor,T1,T2,FA_scale,neg_phase,mask_file)

if ~exist('T1','var') || isempty(T1)
    T1=[2.6,1.2];  %blood, wm.
end

if ~exist('T2','var') || isempty(T2)
    T2=[28.5,23.5]*0.001;  %putaman T2* is 17.6 ms
end


if ~exist('neg_phase','var')
    neg_phase=false;
end

FA_scale=get(params,'flip angle scale');
meas=get_fpattern(params,'protocol');
f_vessel=get_fpattern(params,'vessel mask');
f_wm=get_fpattern(params,'WM mask');

interp_factor=get(params,'interp factor');
bg_size=get(params, 'bg size');
method=get(params,'method');
num2deg=get(params,'num2deg');
mag_file=get_fpattern(params,'mean mag file');
ph_file=get_fpattern(params,'mean phase file');

f_nonzero=get(params,'nonzero fraction');


thk=readsPar(meas,'asSlice[0].dThickness');
thk=thk/10;

TR=readsPar(meas,'alTR[0]');
seg=readsPar(meas,'lSegments');

TR=TR/2/seg/1e6;

TE = readsPar(meas,'alTE[0]');

TE=TE/1e6;

d_blood=1.06;  %g/ml
lambda = 1/0.9; %g/ml; blood to brain water partition;
lambda=lambda/d_blood;  %ml/ml;

dReadoutFOV=readsPar(meas,'asSlice[0].dReadoutFOV');
dPhaseFOV=readsPar(meas,'asSlice[0].dPhaseFOV');

lRO=readsPar(meas,'lBaseResolution');

lPE = readsPar(meas,'lPhaseEncodingLines');

vox_size=[dReadoutFOV,dPhaseFOV]./round([lRO,lPE]*f_nonzero);

VENC = readsPar(meas,'nVelocity');
bipolar= readsPar(meas,'alFree[21]');
if ~isempty(bipolar)
    VENC=VENC/2;
end

global fit_method
%use_fminsearch = false; 
global flow_pattern
%flow_pattern='plug';

FA=readsPar(meas,'adFlipAngleDegree[0]');

fprintf('TE = %3.1f ms; TR = %3.1f ms; VENC = %d cm/s; FA = %d deg; thickness = %2.1f cm\n',TE*1000,TR*1000,VENC,FA,thk);
fprintf('VoxelSize: %3.2f*%3.2f\n',vox_size);

FA=FA*FA_scale;

f_prep=['TR_',num2str(TR),'_T1_',num2str(T1),'_T2_',num2str(T2),'_TE_',num2str(TE),'_thk_',num2str(thk),'_FA_',num2str(FA),'_',flow_pattern,'_sinc_fl_fq_retroz_mb.mat'];
f_prep=strrep(f_prep,' ','_');
for i=1:10
    f_prep=strrep(f_prep,'__','_');
end

d_prep=fullfile(to7t,'pTX_Human','PVS_R21','Prep_PVFlow',f_prep);

if ~exist(d_prep,'file')
    [v_calc,sart_calc,v_calc_pc,sart_calc_pc,s_static]=Prep_PV_V4TOF2Sur_PC(TR,T1,TE,T2,thk,FA,flow_pattern,'sinc','fl_fq_retroZ_mb');
    save(d_prep,'v_calc','sart_calc','v_calc_pc','sart_calc_pc','s_static');
else
    load(d_prep);
end


%%


try
    roi=ri(f_vessel,'','','d');
catch
    roi=ri(f_vessel);
end

try
    m_wm=ri(f_wm,'','','d');
catch
    m_wm=ri(f_wm);
end

% %MID75
% ph_file='Phase_MID75_interp5_mean.mat';
% mag_file='Mag_MID75_interp5_mean.mat';
% roi=ri('mask_vessels_MID75.mat','','','d');
% m_wm=ri('roi_wm_MID75_interp5.mat');
% _us0_75_interp5_mean

if max(roi(:))==1  %not yet clusterized
    roi=clusterize2(roi);
end

nroi=max(roi(:));
pv=zeros(1,nroi);
v=zeros(1,nroi);
rad = zeros(1,nroi);
nv=zeros(1,nroi);
ph_total=zeros(1,nroi);
ph_max=zeros(1,nroi);

center=zeros(2,nroi);
tof=zeros(1,nroi);
ok=zeros(1,nroi);  % is the residual zero?
area=prod(vox_size./interp_factor);
mbg=zeros(2,nroi);

mag=ri_d1(mag_file);
ph=ri_d1(ph_file);

mag0=double(mag);
ph0=double(ph);

if neg_phase
    ph0=-ph0;
end


                         
    do_detrend=true;
for i=1:max(roi(:))
    %%
    
    tic;
    roi2=roi==i;
    nv(i)=sum(roi2(:));
    if nv(i)==0
        continue;
    end
    
    bg_size_i=round(mean(bg_size./vox_size.*interp_factor));
    
     bg=mask_circle(size(roi),bg_size_i,roiCOM(roi==i),1);
     
   % bg=roiCOMNeighbor(roi==i,bg_size_i);
    
    bgroi=bg&roi==0&m_wm>0;
    fprintf('%d bg voxels = %d \n',i,sum(bgroi(:)));
    
    %mbg=bg4medianFilter(mag(:,:,:,1),bg_size_i*3,bgroi);
    %phbg=bg4medianFilter(ph,bg_size_i*3,bgroi);
    
    phbg=0;
    pos=roiCOM(roi2);
    roiRad=get(params,'roi radius (mm)');
    roiRad_i=round(mean(roiRad./vox_size.*interp_factor));  % 1 and 2 should be equal
                 
    cropr=pos(1)-roiRad_i:pos(1)+roiRad_i;       
    cropc=pos(2)-roiRad_i:pos(2)+roiRad_i;
            
    if do_detrend  % do we need to do detrend again?
        if size(ph0,4)>1 % phase images
        ph=tof_bg_rm(ph0(:,:,:,2*i-1:2*i),bgroi,[],[],false);
        mag=tof_bg_rm(mag0(:,:,:,2*i-1:2*i),bgroi);
        else % phase contrast images
         ph=tof_bg_rm(ph0,bgroi,[],[],false);
         mag=tof_bg_rm(mag0,bgroi);
        end
    else
        ph=ph0;
        mag=mag0;
    end
    
    
    
    mbg(:,i)=mean_roi(mag,bgroi>0);
    
    
    
    tof(i)=mean_roi(mag(:,:,:,1),roi2)./mbg(1,i);
    
    ph_total(i)=(mean(ph(roi2))-phbg)*num2deg;
    ph_max(i)=(max(ph(roi2))-phbg)*num2deg;
    
     
    if strcmp(method(1),'3')
                   
        s_static=s_static(1);
        
      %  mbg1=bg4medianFilter(mag(:,:,:,1),bg_size_i*3,bgroi);
      %  mbg2=bg4medianFilter(mag(:,:,:,2),bg_size_i*3,bgroi);
        
        snorm=mean(mbg(:,i))*sart_calc(1)/s_static(1)*lambda;
        
        mag1=mag(:,:,:,1);
        mag2=mag(:,:,:,2);
        
        
        if size(ph,4)==1
           cd=mag2.*exp(1i*(ph-phbg)/180*pi*num2deg) - mag1;% negative phase      
        else
           cd=mag2.*exp(-1i*(ph(:,:,:,2)-phbg)/180*pi*num2deg) - mag1.*exp(-1i*(ph(:,:,:,1)-phbg)/180*pi*num2deg);% negative phase         
        end
        cd=cd/snorm;
        
        z(1)=mean(cd(roi2));
        
        %z(1)=mean(mag2(roi2).*exp(1i*(ph(roi2)-phbg)/180*num2deg*pi)-mag1(roi2));
        
        z(2)=mean(mag1(roi2))/snorm;
        
        %z=z/snorm;
        
        s0=sart_calc(1);
        sart_calc=sart_calc/s0;
        sart_calc_pc=sart_calc_pc/s0;
        s_static=s_static/s0;
        
        if strcmp(method(2),'a')
            [pv(i),v(i),ok(i)]=PV_V4ComplexDiff_TOF(z, v_calc,sart_calc,v_calc_pc,sart_calc_pc,s_static,VENC,lambda,flow_pattern,false);
        elseif strcmp(method(2),'b')
            [pv(i),v(i),ok(i)]=PV_V4ComplexDiff_TOF(z, v_calc,sart_calc,v_calc_pc,sart_calc_pc,s_static,VENC,lambda,flow_pattern,true);
        elseif strcmp(method(2),'c')
            roi3=mask_circle(size(cd),roiRad_i,pos,1);
            data(:,:,i)=cd(cropr,cropc);
            roi3=roi3(cropr,cropc);
            % fitfunc2 = @(x) fitfunc(x(3),x(4),x(1:2),data(:,:,i),roi3(cropr,cropc),FOV,vox_size,vox_size./interp_factor,v_calc_pc,sart_calc_pc,VENC);
            c1=1:size(data,1);
            c2=1:size(data,2);
            
            
            
            
            for ifit =1
                
                if ifit==2
                    icenter=round(center(:,i)./vox_size(1)*interp_factor);
                    c1=(5+icenter(1))+(-2:2);
                    c2=(5+icenter(2))+(-2:2);
                end
                
                mag_save(:,:,i,:)=mag(cropr,cropc,1,:);
                mag_save(:,:,i,:)=setv_roi(mag_save(:,:,i,:),roi3==0,0);
                
                ph_save(:,:,i,:)=ph(cropr,cropc,1,:);
                ph_save(:,:,i,:)=setv_roi(ph_save(:,:,i,:),roi3==0,0);
                
                
                
                FOV = size(data(c1,c2,i)).*vox_size./interp_factor;
                fitfunc2 = @(x) fitfunc(x(3),x(4),x(1:2),data(c1,c2,i),roi3(c1,c2),FOV,vox_size,vox_size./interp_factor,v_calc_pc,sart_calc_pc,VENC,flow_pattern);
                
               skip_fit=0;
               if skip_fit
                   res=[0,0,0,0];
               else
                   
                   if  strcmp(fit_method,'patternsearch')
                      
                   %    rerun= strcmp(output_tmp(i).message,'Maximum number of function evaluations exceeded: increase options.MaxFunctionEvaluations.');
                         
                       options = optimoptions('patternsearch','Display','iter','MaxIter',Inf,'MaxFunctionEvaluations',Inf);
                      rerun=1;
                      if rerun
                     [res,tmp,exitflag(i),output(i)]=patternsearch(fitfunc2,[0,0,0.08,0.6],[],[],[],[],[-0.5,-0.5,0.01,0],[0.5,0.5,0.25,5],[],options);
                     
                      else
                         output(i)=output_tmp(i);
                         exitflag(i)=exitflag_tmp(i);
                         res(4)=v_tmp(i);
                         res(3)=rad_tmp(i);
                         res(1:2)=center_tmp(:,i);
                      end
                   elseif strcmp(fit_method,'fminsearch')
                       options = optimset('Display','iter','MaxIter',1000);
                       res=fminsearch(fitfunc2,[0.5,0.5,0.25,5],options); 
                   else
                       options = optimoptions('lsqnonlin','MaxIter',Inf,'MaxFunctionEvaluations',Inf,'FunctionTolerance',0,'OptimalityTolerance',0);
                       [res,tmp,tmp2,exitflag(i),output(i)]=lsqnonlin(fitfunc2,[0,0,0.08,0.6],[-0.5,-0.5,0.01,0],[0.5,0.5,0.25,5],options);
                   end
               end
               
                if ifit==1
                    center(:,i) = res(1:2);
                end
               
                
            end
            
            pv(i) = pi*res(3)^2/nv(i)/area;
            v(i) = res(4);  % in cm/s
            rad(i)=res(3); %% in mm
            
            % x=[rad,v];
            
            
            [resid,cd_fit(:,:,i)]=fitfunc2(res);
            
            plot_fit=0;
            if plot_fit  %check fit
                disp(res);
                tmp1=data(c1,c2,i);
                tmp1(~roi3)=0;
                tmp2=cd_fit(:,:,i);
                tmp2(~roi3)=0;
                
               figure; compare_cd(tmp1,tmp2);
                resid=sos(tmp1(:)-tmp2(:),1).^2;
                disp(resid);
               return; 
            end
            
            roi_pattern=roi3;
        else
            error('unknown method');
        end
        
    end
    
   fprintf('Time remaining %4.3f s\n',toc*(max(roi(:))-i));
end

    flow=pv.*v*area.*nv*10;  %in mm3/s
    
    flow_uc=ph_total*VENC/180*area.*nv*10;
    flow_max=ph_max*VENC/180*area*10;
    
    dir_name=fileparts(ph_file);
    
    if do_detrend
        
        fname=sprintf('FlowPartialVolume_%s_method%s_oldCD.mat',strtok(filename(f_vessel),'.'),method);
        
        if strcmp(flow_pattern,'plug')
        fname=filename_append(fname,'_plug');
        end
        
        if strcmp(fit_method,'fminsearch')
         fname=filename_append(fname,'_fminsearch');
        elseif strcmp(fit_method,'patternsearch')
            fname=filename_append(fname,'_patternsearch_maxInf');
        else
            fname=filename_append(fname,'_lsqnonlin_OptimFuncTol_0');
        end
        
  
    else
        fname=sprintf('FlowPartialVolume_%s_method%s_no2ndBGdetrend.mat',strtok(filename(f_vessel),'.'),method);
    end
    
    if ~strcmp(method,'3c')
      save(fname,'flow', 'flow_uc', 'pv', 'v', 'ok', 'ph_total', 'tof','nv','area','flow_max','mbg','params','do_detrend');
    else
        save(fname,'params','ph_save','mag_save','flow', 'flow_uc', 'pv', 'v', 'ok', 'ph_total', 'tof','nv','area','flow_max','cd_fit','data','roi_pattern','center','rad','mbg','do_detrend','exitflag','output');
    end
    fprintf('flow (mm3/s) = \n');
    disp([1:length(flow);flow]);
    
    fprintf('mean (sem) flow (mm3/s) = %f(%f)\n',mean(flow),std(flow)/sqrt(length(flow)));
   
            
            
      function res=fr(r,vmean,va,sa,rad,venc,flow_pattern)
    if strcmp(flow_pattern,'laminar')
      vi=v_laminarFlow(r/rad,vmean);
    else
        vi=vmean*ones(size(r));
    end
    si=interp1(va,sa,vi);
    res=si.*(exp(1i*vi/venc*pi));
    res(r>rad)=0;
    
    
        function [res,cd_fit]=fitfunc(rad,vmean,center,cd,roi,FOV,voxSize,voxSize_interp,va,sa,venc,flow_pattern)
            
            fr2=@(r) fr(r,vmean,va,sa,rad,venc,flow_pattern);
            im1=circleImage_CartesianKSpace(fr2,center,FOV,voxSize,voxSize_interp);
            
            fr2=@(r) fr(r,vmean,va,sa,rad,Inf,flow_pattern);
            im0=circleImage_CartesianKSpace(fr2,center,FOV,voxSize,voxSize_interp);
            
            ph=angle(im1.*conj(im0));
       %     cd_fit = abs(im1).*exp(1i*ph)-abs(im0);
            cd_fit=im1-im0;    
        %    cd_fit = abs(im0).*(exp(1i*ph)-1);
            resid=cd(roi>0)-cd_fit(roi>0);
            res=[real(resid(:));imag(resid(:))];
                   global fit_method
                   if strcmp(fit_method,'fminsearch') || strcmp(fit_method,'patternsearch') 
                       
                      res=sum(res.^2,1);                        
                       
                   end
                   
            