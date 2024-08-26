for i=1:10
    fnames_lp{i} = sprintf('rdk%d_dtr_norm_lp0p1+orig',i);
    fnames{i} = sprintf('rdk%d_dtr_norm+orig',i);
end
 
 
cd('/media/SS_FMRI1/MRI_Analysis/Peng/RDK_12_16_2008/afni');


cd /media/SS_FMRI1/MRI_Analysis/Hua/rdk_12_04_08/afni;


cd /media/SS_FMRI1/MRI_Analysis/KELIRIS004_GeorgeSister/rdk/afni;
clear fnames;
clear fnames_lp;
scans = [1,2,3,5,6,7,8];
for i=1:length(scans)
    fnames_lp{i} = sprintf('rdk%d_dtr_norm_lp0p1+orig',scans(i));
    fnames{i} = sprintf('rdk%d_dtr_norm+orig',scans(i));
end

bn= maskcalc('ccmask_lp0p1_p0p05+orig','step(a)',21);
figure;
tsroi_sortave(fnames,'stim/stim_cohall.1D',20,5,4,maskcalc('maskcalc_temp+orig','step(2.5-a)*step(a)'));

plot_labels('V5 && p>0.005','t (TR)','signal change',{'coherence = 25%','50%','100%'},16);

gtext('104 voxels','FontSize',16);


arearoi_sortave(fnames,'stim/stim_cohall.1D',20,5,t_peak,mask,t_base)