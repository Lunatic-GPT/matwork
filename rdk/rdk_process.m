% retinotopy preprocess
flist = {'Ret_ring1','Ret_ring2','Ret_ring3','Ret_wedge1','Ret_wedge2','Ret_wedge3','anatomy'};
acquisition_time_dicom(flist);
afni_preprocess_baylor(flist(1:end-1),'anatomy');

!3dcalc -a s1_dtr_norm+orig'[2..193]' -b s2_dtr_norm+orig'[2..193]' -c s3_dtr_norm+orig'[2..193]' -expr '(a+b+c)/3' -prefix ring_dtr_norm
!3dcalc -a s4_dtr_norm+orig'[2..193]' -b s5_dtr_norm+orig'[2..193]' -c s6_dtr_norm_194brk+orig'[2..193]' -expr '(a+b+c)/3' -prefix wedge_dtr_norm

%rdk preprocess
flist = {'RFLoc','RDK1','RDK2','RDK3','RDK4','RDK5','RDK6','RDK7','RDK8','RDK9','RDK10','anatomy'};
acquisition_time_dicom(flist);
afni_preprocess_baylor(flist(1:end-1),'anatomy');

% align the time courses.
rkd2_gen_seq('../stimulus/*RDKMo*.mat');
p = parameter('sort time series from different trials');
p = add(p,'string','time series','s2_dtr_norm+orig.HEAD,s3_dtr_norm+orig.HEAD,s4_dtr_norm+orig.HEAD,s5_dtr_norm+orig.HEAD,s6_dtr_norm+orig.HEAD,s7_dtr_norm+orig.HEAD,s8_dtr_norm+orig.HEAD,s9_dtr_norm+orig.HEAD,s10_dtr_norm+orig.HEAD,s11_dtr_norm+orig.HEAD');
p = add(p,'string','subbrik range','9..128'); 
p = add(p,'filename','stimulus sequence file','stim_seq.1D');
p = add(p,'int','trial duration (TR)',15);
p = add(p,'int','baseline duration (TR)',0);  
p = add(p,'string','prefix for output files','ts_sort');
tsAlign(p);

%% correlation analysis for the localizer
copyfile /home/zong/rdk2/Sub10_CTR_RDK_102510/Sub10_CTR_RDK_102510/afni/ref_gamma_loc.1D .
correlateAnalysis('s1_dtr_norm+orig','ref_gamma_loc.1D','cc_loc');



