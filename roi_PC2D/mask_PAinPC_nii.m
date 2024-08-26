function mask_PAinPC_nii(fmag,fphase,wm_mask,alpha_2tail,prefix)
% mask_PAinPC_new(fphase,fmag,wm_mask,alpha_2tail)
% calcu a PA mask based on a series of phase and mag images givin by sprintf(fpat,scanIDs(i))
% In the image (e.g. dicom dir), first component is phase, second component mag.
% tSNR calculated to determine threshold for roi determination.
% par should contain a parameter named alpha_2tail 

if ~exist('alpha_2tail','var')
    alpha_2tail=0.1;
end


warning off;

ph=ri(fphase);
mag=ri(fmag);
mag=double(mag);

ph=double(ph)/18000*pi;

nii=load_untouch_niigz(fphase);

voxsize = nii.hdr.dime.pixdim(2:4);

roi=mask_PAinPC(mag,ph,wm_mask,alpha_2tail,voxsize);


nii.img=roi;

save_untouch_niigz(nii,prefix); 
