function ts_sort_bin_gui


rdkStimFile('../RDK_out','stim');
ref_gamma(10,10,1);
conv_stim_ir('refwave.dat','stim/stim_cohall.1D','ref_gamma.1D');   % prepare for glm analysis


for i=1:10
 flist{i} = sprintf('s%d_vdtr_norm+orig',i);
end
nTR_trial = 20;
nTR_stim = 5;
noff1 = 0;


%tssort_bin(flist,'stim_allscans/stim_cohall.1D',nTR_trial,nTR_stim,noff1,'ts_sort'); %K004
tssort_bin(flist,'stim/stim_cohall.1D',nTR_trial,nTR_stim,noff1,'ts_sort_vd');