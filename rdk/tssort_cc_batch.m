for i=1:10
    fnames_lp{i} = sprintf('rdk%d_dtr_norm_lp0p1+orig',i);
    fnames{i} = sprintf('rdk%d_dtr_norm+orig',i);
end
 
 
cd('/media/SS_FMRI1/MRI_Analysis/Peng/RDK_12_16_2008/afni');
tssort_bin(fnames_lp,'stim/stim_cohall.1D',20,5,0,'tssorted_lp0p1');
tssort_bin(fnames,'stim/stim_cohall.1D',20,5,0,'tssorted');


for i=1:4
 correlateAnalysis(['tssorted_lp0p1_coh',num2str(i),'+orig'],'ref_gamma_20.1D',['cc_lp0p1_coh',num2str(i)]);
 correlateAnalysis(['tssorted_coh',num2str(i),'+orig'],'ref_gamma_20.1D',['cc_coh',num2str(i)]);
end

cd /media/SS_FMRI1/MRI_Analysis/Hua/rdk_12_04_08/afni;
tssort_bin(fnames_lp,'stim/stim_cohall.1D',20,5,0,'tssorted_lp0p1');
tssort_bin(fnames,'stim/stim_cohall.1D',20,5,0,'tssorted');


for i=1:4
 correlateAnalysis(['tssorted_lp0p1_coh',num2str(i),'+orig'],'ref_gamma_20.1D',['cc_lp0p1_coh',num2str(i)]);
 correlateAnalysis(['tssorted_coh',num2str(i),'+orig'],'ref_gamma_20.1D',['cc_coh',num2str(i)]);
end

cd /media/SS_FMRI1/MRI_Analysis/KELIRIS004_GeorgeSister/rdk/afni;
clear fnames;
clear fnames_lp;
scans = [1,2,3,5,6,7,8];
for i=1:length(scans)
    fnames_lp{i} = sprintf('rdk%d_dtr_norm_lp0p1+orig',scans(i));
    fnames{i} = sprintf('rdk%d_dtr_norm+orig',scans(i));
    Fourier3d(fnames{i},0.1,fnames_lp{i});  
end
tsave_rdk(fnames_lp,'stim/stim_cohall.1D',20,5,'tsave_lp0p1');
tssort_bin(fnames_lp,'stim/stim_cohall.1D',20,5,0,'tssorted_lp0p1');


tssort_bin(fnames,'stim/stim_cohall.1D',20,5,0,'tssorted');


for i=1:3
 correlateAnalysis(['tssorted_lp0p1_coh',num2str(i),'+orig'],'ref_gamma_20.1D',['cc_lp0p1_coh',num2str(i)]);
% correlateAnalysis(['tssorted_coh',num2str(i),'+orig'],'ref_gamma_20.1D',['cc_coh',num2str(i)]);
end
