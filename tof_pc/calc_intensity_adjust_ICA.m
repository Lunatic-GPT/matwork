function [int_adjust,snorm]=calc_intensity_adjust_ICA(params,vmean,T1,T2)
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

FA_scale=1;%get(params,'flip angle scale');
meas=get_fpattern(params,'protocol');
f_vessel=get_fpattern(params,'vessel mask');
f_wm=get_fpattern(params,'WM mask');

interp_factor=get(params,'interp factor');
mag_file=get_fpattern(params,'mag file');


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

VENC = readsPar(meas,'nVelocity');
bipolar= readsPar(meas,'alFree[21]');
if ~isempty(bipolar)
    VENC=VENC/2;
end

flow_pattern='laminar';

FA=readsPar(meas,'adFlipAngleDegree[0]');

fprintf('TE = %3.1f ms; TR = %3.1f ms; VENC = %d cm/s; FA = %d deg; thickness = %2.1f cm\n',TE*1000,TR*1000,VENC,FA,thk);


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
    roi=ri_d1(f_vessel,'','','d');


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

mag=ri_d1(mag_file);

mag=double(mag);
%roi=clusterize2(roi);

    old_style=0;
    
    if old_style
    m_wm=ri_d1(f_wm);
    bgroi=m_wm>0;
    
    mbg=mean_roi(mag,bgroi>0);
    end
    
for i=1:max(roi(:))
    
    tof(i)=mean_roi(mag(:,:,:,1),roi==i);
    
    sart=interp1(v_calc,sart_calc,abs(vmean(i)));
    
    if old_style
        snorm_tmp=mean(mbg(:,i))*sart_calc(1)/s_static(1)*lambda; % should be
        measured=tof(i)/snorm_tmp;%
        
        theory=sart/sart_calc(1);    
        int_adjust(i)= theory/measured;  % no used any more
    else
        int_adjust(i)=1;  % not use, set to any value
    end
        
    
    snorm(i)=tof(i)*sart_calc(1)/sart;
   %     snorm=mean(mbg(:,i))*sart_calc(1)/s_static(1)*lambda/int_adjust;
        
     
end

            




            
            
            
            