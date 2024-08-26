function Flow_PartialVolume4Pulsatility(params,T1,T2,neg_phase)
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
mag_file=get_fpattern(params,'mag file');
ph_file=get_fpattern(params,'phase file');

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

vox_size=[dReadoutFOV,dPhaseFOV]./[lRO,lPE];

VENC = readsPar(meas,'nVelocity');
bipolar= readsPar(meas,'alFree[21]');
if ~isempty(bipolar)
    VENC=VENC/2;
end
nvenc=length(VENC)+1;

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

area=prod(vox_size./interp_factor);


mag=ri_d1(mag_file);
ph=ri_d1(ph_file);

mag0=double(mag);
ph0=double(ph);
nt=size(ph0,4);


nroi=max(roi(:));
v=zeros(nroi,nt);

mbg=zeros(nvenc,nroi);


if neg_phase
    ph0=-ph0;
end

%roi=clusterize2(roi);

dir_name=fileparts(ph_file);
if ~isempty(strfind(f_vessel,'tof_vmask'))
    fname=fullfile(dir_name,sprintf('Results_Flow_PartialVolume_method%s_rmax0_15mm_tof_vmask.mat',method));
elseif ~isempty(strfind(f_vessel,'Frangi_vmask'))
    fname=fullfile(dir_name,sprintf('Results_Flow_PartialVolume_method%s_rmax0_15mm_Frangi_vmask.mat',method));
else
    fname=fullfile(dir_name,sprintf('Results_Flow_PartialVolume_method%s_rmax0_15mm.mat',method));
end

x0=load(fname,'rad','pv','v','center');
if max(roi(:))~=length(x0.v)
    error('max(roi(:))~=length(x0.v)');
    
end
for i=1:max(roi(:))
    %%
    tic;
    roi2=roi==i;
    nv(i)=sum(roi2(:));
    if nv(i)==0
        continue;
    end
    
    bg_size_i=round(mean(bg_size./vox_size.*interp_factor));
    
    bg=roiCOMNeighbor(roi==i,bg_size_i);
    
    bgroi=bg&roi==0&m_wm>0;
    fprintf('%d bg voxels = %d \n',i,sum(bgroi(:)));
    
    
    phbg=0;
    ph=tof_bg_rm(ph0,bgroi,[],[],false,true);
    
    mag0=reshape(mag0,[size(mag0(:,:,1,1)),nvenc,nt]);
    mag=0*mag0;
    for j=1:nvenc
        mag(:,:,j,:)=tof_bg_rm(mag0(:,:,j,:),bgroi,[],[],true,true);
        
    end
    
    mbg(:,i)=mean_roi(mean(mag,4),bgroi>0);
    
    
    
    if strcmp(method(1),'3')
        s_static=s_static(1);
        
        %  mbg1=bg4medianFilter(mag(:,:,:,1),bg_size_i*3,bgroi);
        %  mbg2=bg4medianFilter(mag(:,:,:,2),bg_size_i*3,bgroi);
        
        snorm=mean(mbg(:,i))*sart_calc(1)/s_static(1)*lambda;
        
        
        for it =1:nt
            mag1=mag(:,:,1,it);
            mag2=mag(:,:,2,it);
            
            cd=mag2.*exp(1i*(ph(:,:,:,it)-phbg)/180*pi*num2deg)/snorm - mag1/snorm;% negative phase
            z(1)=mean(cd(roi2));
            
            %z(1)=mean(mag2(roi2).*exp(1i*(ph(roi2)-phbg)/180*num2deg*pi)-mag1(roi2));
            
            z(2)=mean(mag1(roi2))/snorm;
            
            %z=z/snorm;
            
            s0=sart_calc(1);
            sart_calc=sart_calc/s0;
            sart_calc_pc=sart_calc_pc/s0;
            s_static=s_static/s0;
            
            if strcmp(method(2),'c')
                
                roiRad=get(params,'roi radius (mm)');
                roiRad_i=round(mean(roiRad./vox_size.*interp_factor));  % 1 and 2 should be equal
                
                pos=roiCOM(roi2);
                roi3=mask_circle(size(cd),roiRad_i,pos,1);
                
                cropr=pos(1)-roiRad_i:pos(1)+roiRad_i;
                
                cropc=pos(2)-roiRad_i:pos(2)+roiRad_i;
                
                data=cd(cropr,cropc);
                
                FOV = size(data).*vox_size./interp_factor;
                
                fitfunc2 = @(x) fitfunc(x0.rad(i),x,x0.center(1:2,i),data,ones(size(data)),FOV,vox_size,vox_size./interp_factor,v_calc_pc,sart_calc_pc,VENC);
                
                v(i,it)=lsqnonlin(fitfunc2,1.3538,0,VENC); % in cm/s
                
                
                [resid,cd_fit(:,:,i,it)]=fitfunc2(v(i,it));
                
            else
                error('unknown method');
            end
        end
    end
fprintf('time remaining %f s\n',toc*(max(roi(:))-i));
end

save(filename_append(fname,'_pulsatility'),'v');


%%


function res=fr(r,vmean,va,sa,rad,venc)

vi=v_laminarFlow(r/rad,vmean);
si=interp1(va,sa,vi);
res=si.*(exp(1i*vi/venc*pi)-1);
res(r>rad)=0;


function [res,cd_fit]=fitfunc(rad,vmean,center,cd,roi,FOV,voxSize,voxSize_interp,va,sa,venc)

fr2=@(r) fr(r,vmean,va,sa,rad,venc);
cd_fit=circleImage_CartesianKSpace(fr2,center,FOV,voxSize,voxSize_interp);

resid=cd(roi>0)-cd_fit(roi>0);
res=[real(resid(:));imag(resid(:))];

