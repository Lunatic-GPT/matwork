function Flow_PartialVolume_ICA(params,T1,T2,neg_phase)
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


int_adjust=get(params,'intensity adjust');  % factor for adjusting intensity inhomogeneity; intensity ratio: WM:ICA
snorm_external=get(params,'snorm');  % factor for adjusting intensity inhomogeneity; intensity ratio: WM:ICA

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
num2deg=get(params,'num2deg');
mag_file=get_fpattern(params,'mean mag file');
ph_file=get_fpattern(params,'mean phase file');

f_nonzero=get(params,'nonzero fraction');

method = get(params,'method');
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

flow_pattern='laminar';

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

% %% test k-space undersample instead of convolution.
% venc=90;
% kmax=pi/4.5;
% kx=linspace(-kmax,kmax,65);
% kx=kx(1:end-1);
% 
% ky=linspace(-kmax,kmax,65);
% ky=ky(1:end-1);
% 
% kx=repmat(kx,[length(ky),1]);
% ky=repmat(ky(:),[1,length(kx)]);
% 
% fr=@(x) Signal_laminarFlow(x,20.8,v_calc_pc,sart_calc_pc,2.16,venc);
% kdata=kspaceData_circleImage(fr,false,20,2.16,kx,ky);
% kdata=zpad(kdata,640,640);
% im2=ifft2c(kdata);
% 
% fr=@(x) Signal_laminarFlow(x,20.8,v_calc_pc,sart_calc_pc,2.16,Inf);
% kdata=kspaceData_circleImage(fr,false,20,2.16,kx,ky);
% kdata=zpad(kdata,640,640);
% im1=ifft2c(kdata);
% 
% cd=im2-im1;
% cd2=cd(321-19:320+19,321-19:320+19);
% im=cat(3,real(cd2),imag(cd2));
% 
% figure;imshow4(im,[],[1,2]);
% 
% figure;imshow(real(cd2),[]);
% figure;imshow(imag(cd2),[]);




%%

%%
try
    roi=ri(f_vessel,'','','d');
catch
    roi=ri(f_vessel);
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

%method 
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

%roi=clusterize2(roi);

for i=1:max(roi(:))
    %%
    
    tic;
    roi2=roi==i;
    nv(i)=sum(roi2(:));
    if nv(i)==0
        continue;
    end
    
    
    
    phbg=0;
ph=ph0;
mag=mag0;


  
    
    ph_total(i)=(mean(ph(roi2))-phbg)*num2deg;
    ph_max(i)=(max(ph(roi2))-phbg)*num2deg;
    
        
        s_static=s_static(1);
        
        
        mag1=mag(:,:,:,1);
        mag2=mag(:,:,:,2);
        ph1=ph(:,:,:,1);
        ph2=ph(:,:,:,2);
        
        
        if isempty(snorm_external)  %not used any more
            
            m_wm=ri_d1(f_wm,'','','d');
            bgroi=m_wm>0;
            mbg(:,i)=mean_roi(mag,bgroi>0);
            snorm=mean(mbg(:,i))*sart_calc(1)/s_static(1)*lambda/int_adjust(i);  %snorm = mag*sart(0)/sart(v);
        else
            snorm=snorm_external(i);
        end
        
      %  cd=mag2.*exp(1i*(ph-phbg)/180*pi*num2deg)/snorm - mag1/snorm;%
      %  this will cause problem because mag1 is negative for the first
      %  side lobe, therefore use mag2*exp(1i*ph2)-mag1*exp(1i*ph1)
      
      cd=mag2.*exp(-1i*ph2/180*pi*num2deg)/snorm-mag1.*exp(-1i*ph1/180*pi*num2deg)/snorm;
        
    %tof(i)=mean_roi(mag(:,:,:,1),roi2)./mbg(1,i);
      %  cd =do_undersample_interp(cd,0.1,10);
%         
%         z(1)=mean(cd(roi2));
%         
%         %z(1)=mean(mag2(roi2).*exp(1i*(ph(roi2)-phbg)/180*num2deg*pi)-mag1(roi2));
%         
%         z(2)=mean(mag1(roi2))/snorm;
        
        %z=z/snorm;
        
        s0=sart_calc(1);
        sart_calc=sart_calc/s0;
        sart_calc_pc=sart_calc_pc/s0;
        s_static=s_static/s0;
        
            
            roiRad=get(params,'roi radius (mm)');
            roiRad_i=round(mean(roiRad./vox_size.*interp_factor));  % 1 and 2 should be equal
            
            pos=roiCOM(roi2);
            roi3=mask_circle(size(cd),roiRad_i,pos,1);
            
            cropr=pos(1)-roiRad_i:pos(1)+roiRad_i;
            
            cropc=pos(2)-roiRad_i:pos(2)+roiRad_i;
            
            data(:,:,i)=cd(cropr,cropc);
            
            mag_crop(:,:,:,:,i)=mag(cropr,cropc,:,:);
              
            ph_crop(:,:,:,:,i) = (ph(cropr,cropc,:,:)-phbg)/180*pi*num2deg;
            
            FOV = size(data(:,:,i)).*vox_size./interp_factor;
           % fitfunc2 = @(x) fitfunc(x(3),x(4),x(1:2),data(:,:,i),roi3(cropr,cropc),FOV,vox_size,vox_size./interp_factor,v_calc_pc,sart_calc_pc,VENC);
             fitfunc2 = @(x) fitfunc(x(3),x(4),x(1:2),data(:,:,i),ones(size(data(:,:,1))),FOV,vox_size,vox_size./interp_factor,v_calc_pc,sart_calc_pc,VENC);
            
          %  res=lsqnonlin(fitfunc2,[0,0,0.0824,1.3538],[-0.5,-0.5,0,0],[0.5,0.5,1,8]);
            res(i,:)=lsqnonlin(fitfunc2,[0,0,2.0,17],[-6,-6,0.4,0],[6,6,6,90]);
            
            %   res(i,:)=lsqnonlin(fitfunc2,[0,0,2.0,17],[-6,-6,2.0,0],[6,6,2.0,90]);
          
            pv(i) = pi*res(3)^2/nv(i)/area;
            v(i) = res(i,4);  % in cm/s
            rad(i)=res(i,3); %% in mm
            center(:,i) = res(i,1:2);
            % x=[rad,v];
            
            
            [resid,cd_fit(:,:,i)]=fitfunc2(res(i,:));
            
            roi_pattern=roi3(cropr,cropc);
             
            
           
        
   fprintf('Time remaining %4.3f s\n',toc*(max(roi(:))-i));
end


    flow=pi*rad.^2.*v*10;  %in mm3/s
    
    flow_uc=ph_total*VENC/180*area.*nv*10;
    flow_max=ph_max*VENC/180*area*10;
    
    dir_name=fileparts(ph_file);
    
    fname=fullfile(dir_name,sprintf('Results_Flow_PartialVolume_method%s_new_roi2.mat',method));
    if ~strcmp(method,'3c')
      save(fname,'flow', 'flow_uc', 'pv', 'v', 'ok', 'ph_total', 'nv','area','flow_max','mbg','res','params','cd_fit','data');
    else
        save(fname,'mag_crop','ph_crop','flow', 'flow_uc', 'pv', 'v', 'ok', 'ph_total', 'nv','area','flow_max','cd_fit','data','roi_pattern','center','rad','mbg','res','params');
    end
    fprintf('flow (cm3/s); rad (mm); v (cm/s) = \n');
    disp([1:length(flow);flow/1000;rad;v]);
    
    fprintf('mean (sem) flow (cm3/s) = %f(%f)\n',mean(flow)/1000,std(flow)/1000/sqrt(length(flow)));
    %%
    
%% calculate surrogate CD of undersampled images.

    function res=fr(r,vmean,va,sa,rad,venc)
    
    vi=v_laminarFlow(r/rad,vmean);
    si=interp1(va,sa,vi);
    res=si.*(exp(1i*vi/venc*pi));
    res(r>rad)=0;
    
    
        function [res,cd_fit]=fitfunc(rad,vmean,center,cd,roi,FOV,voxSize,voxSize_interp,va,sa,venc)
            
            fr2=@(r) fr(r,vmean,va,sa,rad,venc);
            im1=circleImage_CartesianKSpace(fr2,center,FOV,voxSize,voxSize_interp);
            
            fr2=@(r) fr(r,vmean,va,sa,rad,Inf);
            im0=circleImage_CartesianKSpace(fr2,center,FOV,voxSize,voxSize_interp);
            
           % ph=angle(im1.*conj(im0));
          %  cd_fit = abs(im1).*exp(1i*ph)-abs(im0);
            cd_fit=im1-im0;
        %    cd_fit = abs(im0).*(exp(1i*ph)-1);
            resid=cd(roi>0)-cd_fit(roi>0);
            res=[real(resid(:));imag(resid(:))];
            
            
            




            
            
            
            