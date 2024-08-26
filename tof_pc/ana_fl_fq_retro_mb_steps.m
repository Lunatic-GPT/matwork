%1. first recon (oversampled data)

fsb='meas_MID71_fl_tof_3slice_105V_FID14273.dat';
fmb='meas_MID69_fl_fq_retroZ_mb_Sinc_FID14271.dat';

recon_fl_fq_mb_grappa(fmb,fsb);


%% 2
prefix=strtok(fmb,'.');

fl_fq_pc(['PC_',prefix,'_mean.mat'],0.04,['Mag_',prefix,'_mean.mat']);

fl_fq_pc(['PC_',prefix,'.mat'],0.04,['Mag_',prefix,'_mean.mat']);

%%
% 3. draw WM ROI with img;


%% 6: detect PVS
pc=['PC_',prefix,'_mean_detrend.mat'];
prefixm='roi_wm_sinc_mb_slice123';
detectPVSCircleInMask(pc,[prefixm,'.mat'],10,3,true);

%%
m=ri([prefixm,'_detectedPVS.mat']);
m_wm=ri([prefixm,'.mat']);
pc=ri(['PC_',prefix,'_detrend_samebg.mat']);
for i=1:size(m,3)
    tmp=clusterize2(m(:,:,i));
    disp([i,max(tmp(:))]);
    
    venc=4;
use_peak=false;
plot_retro_time_course(pc(:,:,i,:),m(:,:,i,:),m_wm(:,:,i,:),venc,use_peak);


end

%%
%%7 show the results



