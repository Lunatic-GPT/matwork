% retinotopy preprocess
%flist = {'Ret_ring1','Ret_ring2','Ret_ring3','Ret_wedge1','Ret_wedge2','Ret_wedge3','anatomy'};
flist = {'ring1','ring2','ring3','wedge1','wedge2','wedge3','anatomy'};
acquisition_time_dicom(flist);
afni_preprocess_baylor(flist(1:end-1),'anatomy');

!3dcalc -a s1_vr+orig'[2..193]' -b s2_vr+orig'[2..193]' -c s3_vr+orig'[2..193]' -expr '(a+b+c)/3' -prefix ring_vr
!3dcalc -a s4_vr+orig'[2..193]' -b s5_vr+orig'[2..193]' -c s6_vr+orig'[2..193]' -expr '(a+b+c)/3' -prefix wedge_vr
!3dDetrend -polort 2 -prefix wedge_vr_dtr wedge_vr+orig
!3dDetrend -polort 2 -prefix ring_vr_dtr ring_vr+orig
retino_wdgrng('wedge_vr_dtr+orig','ring_vr_dtr+orig','retino',32);
%rdk preprocess

flist ={'RDK1','RDK2','RDK3','RDK4','RDK5','RDK6','RDK7','RDK8','RDK9','RDK10','anatomy'};
%flist = {'RFLoc','RDK1','RDK2','RDK3','RDK4','RDK5','RDK6','RDK7','RDK8','RDK9','RDK10','anatomy'};
%flist = {'RFLoc','RDK1','RDK2','RDK3','RDK4','RDK5','RDK6','RDK7','RDK8','anatomy'};
%flist = {'RFLoc0','RDK1','RDK2','RDK3','RDK4','RDK5','RDK6','RDK7','RDK8','RDK9','RDK10','anatomy'};
acquisition_time_dicom(flist);
afni_preprocess_baylor(flist(1:end-1),'anatomy');
afni_preprocess_baylor(flist(1:end-1),'');
!3dWarp -card2oblique s1_mean+orig -prefix anat_wp+orig anat+orig
!@auto_tlrc -base TT_N27+tlrc -input anat_wp+orig

% align the time courses.
rdk2_gen_seq('../stimulus/*RDKMo*.mat');
p = parameter('sort time series from different trials');
%p = add(p,'string','time series','s2_dtr_norm+orig.HEAD,s3_dtr_norm+orig.HEAD,s4_dtr_norm+orig.HEAD,s5_dtr_norm+orig.HEAD,s6_dtr_norm+orig.HEAD,s7_dtr_norm+orig.HEAD,s8_dtr_norm+orig.HEAD,s9_dtr_norm+orig.HEAD,s10_dtr_norm+orig.HEAD,s11_dtr_norm+orig.HEAD');
p = add(p,'string','time series','s2_dtr_norm+orig.HEAD,s3_dtr_norm+orig.HEAD,s4_dtr_norm+orig.HEAD,s5_dtr_norm+orig.HEAD,s6_dtr_norm+orig.HEAD,s7_dtr_norm+orig.HEAD,s8_dtr_norm+orig.HEAD,s9_dtr_norm+orig.HEAD');
p = add(p,'string','subbrik range','9..128'); 
p = add(p,'filename','stimulus sequence file','stim_seq.1D');
p = add(p,'int','trial duration (TR)',15);
p = add(p,'int','baseline duration (TR)',0);  
p = add(p,'string','prefix for output files','ts_sort');
tsAlign(p);
ts_sort2ave('ts_sort_coh1+orig',20,'ts_mean_p1');
ts_sort2ave('ts_sort_coh2+orig',20,'ts_mean_p2');
ts_sort2ave('ts_sort_coh3+orig',20,'ts_mean_p3');
ts_sort2ave('ts_sort_coh4+orig',20,'ts_mean_p4');


ts_sort2ave('ts_sort_p1+orig',15,'ts_mean_p1');
ts_sort2ave('ts_sort_p2+orig',15,'ts_mean_p2');
ts_sort2ave('ts_sort_p3+orig',15,'ts_mean_p3');
ts_sort2ave('ts_sort_p4+orig',15,'ts_mean_p4');
%% correlation analysis for the localizer
!cp /home/zong/rdk2/S10/rdk/afni/ref_gamma_?v.1D .
correlateAnalysis('s1_dtr_norm+orig','ref_gamma_lv.1D','cc_loc_lv');
correlateAnalysis('s1_dtr_norm+orig','ref_gamma_rv.1D','cc_loc_rv');

%% glm analysis
seq2xmat_gamma('stim_seq.1D',9,0,6,9,2,'xmat_gamma.1D');
!cp /home/zong/rdk2/S10/rdk/afni/glm_afni_bg.m .

%% mask generation.
maskcalc('glm_gamma_REML+orig[0]','glm_gamma_REML+orig[7]','(step(a-4.646)+step(-4.646-a))*step(b)',7,'mask_glm');
maskcalc('glm_gamma_REML+orig[0]','step(a-4.646)+step(-4.646-a)',7,'mask_glm');

maskcalc('glm_REML+orig[0]','glm_REML+orig[7]','(step(a-4.619)+step(-4.619-a))*step(b)',7,'mask_glm');
maskcalc('glm_REML+orig[0]','step(a-4.619)+step(-4.619-a)',7,'mask_glm');

%  clean intermediate files

delete_files('s?_tsh+orig.*','.');
delete_files('s??_tsh+orig.*','.');
delete_files('s?_vr_dtr+orig.*','.');
delete_files('s??_vr_dtr+orig.*','.');
delete_files('s?_blur+orig.*','.');
delete_files('s??_blur+orig.*','.');

[nv,ts]=nv_roi('ts_sort_p1+orig','mask_rV5+orig');
plot_sorted_ts(ts,15);




