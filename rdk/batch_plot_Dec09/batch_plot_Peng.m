function batch_plot

% ts Plots
cd /media/Human_fMRI_backup/data/MRI/Peng/RDK_12_16_2008/afni;
root = '/media/Human_fMRI_backup/data/MRI/';

p = parameter('Sort and plot average time series in roi');
p = add(p,'filename','stimulus file',fullfile(root,'Peng/RDK_12_16_2008/afni/stim/stim_cohall.1D'));
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

 %left MT
 mstr = 'mask_rdkR+orig;lMST_thr3p292+orig;equals(a,3)*iszero(b)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'lMT_notMST');

%left MST
mstr = 'lMST_thr3p292+orig;equals(a,2)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'lMST');

%right MT
 mstr = 'mask_rdkL+orig;rMST_thr3p292+orig;equals(a,1)*iszero(b)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'rMT_notMST');

%right MST
mstr = 'rMST_thr3p292+orig;equals(a,1)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'rMST');

%left MT
 mstr = 'mask_rdkR+orig;lMST_thr3p292+orig;mask_rdkL+orig;rMST_thr3p292+orig;equals(a,3)*iszero(b)+equals(c,1)*iszero(d)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'MT_notMST');

%MST
mstr = 'rMST_thr3p292+orig;lMST_thr3p292+orig;equals(a,1)+equals(b,2)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'MST');

% hMT+ from  glm
mstr = 'mask_REML_cohall+orig;equals(a,1)+equals(a,2)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'hMT_glm');
 
% V3a
mstr = 'mask_REML_cohall+orig;equals(a,3)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'V3a');

% ST
mstr = 'mask_REML_cohall+orig;equals(a,4)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'SP');
% rSP
mstr = 'mask_REML_cohall+orig;equals(a,5)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'precentral');
% middle occiptial gyrus
mstr = 'mask_REML_cohall+orig;equals(a,6)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'lST');
mstr = 'mask_REML_cohall+orig;equals(a,7)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'rST');
mstr = 'mask_REML_cohall+orig;equals(a,6)+equals(a,7)';
run_tsroi_arearoi_sortave(p_a,p_ts,mstr,'ST');


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




