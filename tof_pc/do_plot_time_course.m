
mid='MID151';
vox_size=0.3/2.5; 

vessel_mask=sprintf('mask_vessels_%s.mat',mid);
ph_file=sprintf('Phase_%s.mat',mid);
wm_mask=sprintf('mask_wm_%s.mat',mid);
mag_file=sprintf('Mag_%s.mat',mid);
bg_size=round(2/vox_size);

% 
% vessel_mask=ri(vessel_mask,'','','d');
% wm_mask=ri(wm_mask,'','','d');
% ph_file=ri(ph_file);
% 
% 
% vessel_mask=vessel_mask(1:3:end,1:3:end,:,:);
% ph_file=ph_file(1:3:end,1:3:end,:,:);
% wm_mask=wm_mask(1:3:end,1:3:end,:,:);
% bg_size=11;
%ph_file=ri(ph_file);
%ph_file=ph_file/100;

%vessel_mask='mask_surface_vessels.mat';
%mag='Mag_MID75_interp5.mat';  % for SNR estimate
%wm_mask='roi_wm_MID75_interp5.mat';

VENC=4;
use_peak=true;


%[vall,vall_samebg,vall_nobg,vall_bg]=retro_time_course(ph_file,vessel_mask,wm_mask,VENC,use_peak,bg_size,0.01,33*3);
[vall,vall_samebg,vall_nobg,vall_bg]=retro_time_course(ph_file,vessel_mask,wm_mask,VENC,use_peak,bg_size,0.01,[],[-1,1]);


%% estimate velocity measurement error

esl=2;
try
dwm=ri(wm_mask,'','','d');
catch
dwm=ri(wm_mask);
    
end
dph=ri(ph_file);
dmag=ri(mag_file);
VENC=4;
%estimate from magnitude image
val=mean_roi(mean(double(dmag(:,:,esl,1:2:end)),4),dwm(:,:,esl));
eval=mean_roi(std(double(dmag(:,:,esl,1:2:end)),[],4),dwm(:,:,esl));

snr=val/eval*sqrt(12);

%estimate directly from the phase image
eph=mean_roi(std(double(dph(:,:,esl,:)),[],4),dwm(:,:,esl));
ev=eph*VENC/180/100;

disp([VENC/snr/pi*sqrt(2),ev/sqrt(12)]);

%% plot
vsel=1:size(vall,1);
%vsel([17,21])=[];
nshift=5;

h1=figure;
h2=figure;
%vall_136b=circshift(vall_136,[0,1]);
sl=2;
plot_retro_time_course(vall{sl},vall_samebg{sl},vall_nobg{sl},h1,h2,nshift,ev);

    
    
    