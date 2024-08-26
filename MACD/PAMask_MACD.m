function PAMask_MACD(fmag,fph,pro,wm_mask)
% fmag: mag image (not averaget); averaged image should have name _mean appended
% fph: phase image (not averaged);
% scan: scan indices
% wm_mask: the white matter mask
roi_name=mask_PAinPC_new(fmag,fph,wm_mask,0.1);

%% do the fit for individual data and individual roi
par.neg_phase=false;
par.pro=pro;
run_fit_MACD_nii(filename_append(fmag,'_mean',1),filename_append(fph,'_mean',1),[roi_name,'.nii.gz'],par);

