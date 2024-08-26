for i=1:10
    fnames_lp{i} = sprintf('rdk%d_dtr_norm_lp0p1+orig',i);
    fnames{i} = sprintf('rdk%d_dtr_norm+orig',i);
end
 
cd /media/SS_FMRI1/MRI_Analysis/Zong/rdk_12_01_08/afni;
%glm_afni(fnames_lp,'xmat.1D',[],'brainmask+orig','rdk_glm_lp');


glm_afni(fnames,'xmat.1D',[],'brainmask+orig','rdk_glm');


cd('/media/SS_FMRI1/MRI_Analysis/Peng/RDK_12_16_2008/afni');
glm_afni(fnames_lp,'xmat_mc.1D',[],'brainmask+orig','rdk_glm_lp_mc');
glm_afni(fnames,'xmat_mc.1D',[],'brainmask+orig','rdk_glm_mc');


cd /media/SS_FMRI1/MRI_Analysis/Hua/rdk_12_04_08/afni;
glm_afni(fnames_lp,'xmat.1D',[],'brainmask+orig','rdk_glm_lp');
glm_afni(fnames,'xmat.1D',[],'brainmask+orig','rdk_glm');


cd /media/SS_FMRI1/MRI_Analysis/KELIRIS004_GeorgeSister/rdk/afni;
clear fnames;
clear fnames_lp;
scans = [1,2,3,5,6,7,8];
for i=1:length(scans)
    fnames_lp{i} = sprintf('rdk%d_dtr_norm_lp0p1+orig',scans(i));
    fnames{i} = sprintf('rdk%d_dtr_norm+orig',scans(i));
end
glm_afni(fnames_lp,'xmat.1D',[],'brainmask+orig','rdk_glm_lp');
glm_afni(fnames,'xmat.1D',[],'brainmask+orig','rdk_glm');
