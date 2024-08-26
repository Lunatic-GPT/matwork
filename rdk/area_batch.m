
for i=1:10
    fnames{i} = sprintf('rdk%d_dtr_norm_lp0p1+orig',i);
end
 coh =[0.125,0.25,0.5,1];
 
cd('/media/SS_FMRI1/MRI_Analysis/Peng/RDK_12_16_2008/afni');
area_sortave(fnames,'stim/stim_cohall.1D',20,5,1:9,16:20,'area_cohall');
 
 
cd /media/SS_FMRI1/MRI_Analysis/Hua/rdk_12_04_08/afni;
area_sortave(fnames,'stim/stim_cohall.1D',20,5,1:9,16:20,'area_trials');

cd /media/SS_FMRI1/MRI_Analysis/Zong/rdk_12_01_08/afni;
area_sortave(fnames,'stim/stim_cohall.1D',20,5,1:9,16:20,'area_trials');

cd /media/SS_FMRI1/MRI_Analysis/KELIRIS004_GeorgeSister/rdk/afni;
clear fnames;
scans = [1,2,3,5,6,7,8];
coh = [0.25,0.5,1];
for i=1:length(scans)
    fnames{i} = sprintf('rdk%d_dtr_norm_lp0p1+orig',scans(i));
end
area_sortave(fnames,'stim/stim_cohall.1D',20,5,1:9,16:20,'area_trials');

cd /media/SS_FMRI1/MRI_Analysis/KELIRIS005/rdk/afni;
clear fnames;
scans = [1,4,5];
for i=1:length(scans)
    fnames{i} = sprintf('rdk%d_dtr_norm_lp0p1+orig',scans(i));
end
area_sortave(fnames,'stim/stim_cohall.1D',20,5,1:9,16:20,'area_trials');


