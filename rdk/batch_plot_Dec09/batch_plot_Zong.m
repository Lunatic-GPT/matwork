function batch_plot_Zong

%root = '/media/SS_FMRI1/MRI_Analysis/';
cd /media/Human_fMRI_backup/data/MRI/Zong/rdk_12_01_08/afni/;
root = '/media/Human_fMRI_backup/data/MRI/';

p = parameter('Sort and plot average time series in roi');
p = add(p,'filename','stimulus file',fullfile(root,'Zong/rdk_12_01_08/afni/stim/stim_cohall.1D'));
p = add(p,'int','# TR/trial',20);
p = add(p,'int','# TR for stimulus on',5); 
%p = add(p,'int','points for area',5);
p = add(p,'int','# TR before stimulus on',4);
p = add(p,'string','roi mask','maskcalc_temp+orig;equals(a,1)');
p = add(p,'string','save file name','ts.mat');
 
str ='rdk1_dtr_norm+orig';
for i=2:10
     str = sprintf('%s,rdk%d_dtr_norm+orig',str,i);
 end
 p = add(p,'string','time series names',str);

p_a = [];
p_ts = p;

 [mr,nr]=mask_str('brainmask_sroi2v_rh+orig');
 [ml,nl]=mask_str('brainmask_sroi2v_lh+orig');
 
% hMT+ from  glm
mstr = 'mask_glm_REML+orig;equals(a,1)+equals(a,3)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'hMT_glm');
% ST
mstr = 'mask_glm_REML+orig;equals(a,7)+equals(a,2)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'ST');
% rSP
mstr = 'mask_glm_REML+orig;equals(a,4)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'rSP');
% middle occiptial gyrus
mstr = 'mask_glm_REML+orig;equals(a,6)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'lMOG');
mstr = 'mask_glm_REML+orig;equals(a,5)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'rMOG');
mstr = 'mask_glm_REML+orig;equals(a,6)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'lMOG');
mstr = 'mask_glm_REML+orig;equals(a,6)+equals(a,5)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'MOG');

 %V1
mstr= sprintf('%s;%s;%s',mr{1},ml{1},'or(a,b)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V1');

mstr= sprintf('%s;%s;mask_ecc+orig;%s',mr{1},ml{1},'or(a,b)*equals(c,1)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V1_small_ecc');

mstr= sprintf('%s;%s;mask_ecc+orig;%s',mr{1},ml{1},'or(a,b)*equals(c,2)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V1_large_ecc');

%V2
mstr= sprintf('%s;%s;%s;%s;%s',mr{2:3},ml{2:3},'or(a,b,c,d)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V2');
mstr= sprintf('%s;%s;%s;%s;mask_ecc+orig;%s',mr{2:3},ml{2:3},'or(a,b,c,d)*equals(e,1)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V2_small_ecc');
mstr= sprintf('%s;%s;%s;%s;mask_ecc+orig;%s',mr{2:3},ml{2:3},'or(a,b,c,d)*equals(e,2)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V2_large_ecc');
%V3
mstr= sprintf('%s;%s;%s;%s;%s',mr{4:5},ml{4:5},'or(a,b,c,d)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V3');
mstr= sprintf('%s;%s;%s;%s;mask_ecc+orig;%s',mr{4:5},ml{4:5},'or(a,b,c,d)*equals(e,1)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V3_small_ecc');
mstr= sprintf('%s;%s;%s;%s;mask_ecc+orig;%s',mr{4:5},ml{4:5},'or(a,b,c,d)*equals(e,2)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V3_large_ecc');
%V4
mstr= sprintf('%s;%s;%s',mr{6},ml{6},'or(a,b)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V4');
mstr= sprintf('%s;%s;mask_ecc+orig;%s',mr{6},ml{6},'or(a,b)*equals(c,1)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V4_small_ecc');
mstr= sprintf('%s;%s;mask_ecc+orig;%s',mr{6},ml{6},'or(a,b)*equals(c,2)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V4_large_ecc');
%MTLoc
mstr= sprintf('%s;%s;%s',mr{7},ml{7},'or(a,b)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'MTLoc');

mstr= sprintf('%s;%s;mask_MTLoc_ap+orig;%s',mr{7},ml{7},'and(or(a,b),equals(c,1))');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'MTLoc_a');
mstr= sprintf('%s;%s;mask_MTLoc_ap+orig;%s',mr{7},ml{7},'and(or(a,b),equals(c,2))');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'MTLoc_p');


%out rois
mstr= sprintf('%s;%s;%s',mr{8},ml{8},'or(a,b)');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'outside_rois');
%brainmask
mstr= sprintf('brainmask+orig;a');
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'brainmask');



