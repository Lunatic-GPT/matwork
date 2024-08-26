function undersample_image(fname,factor)

nii=load_untouch_niigz(fname);

nii.img=nii.img(1:factor(1):end,1:factor(2):end,1:factor(3):end,:);

nii.hdr.dime.dim(2:4)=nii.hdr.dime.dim(2:4)./factor;
nii.hdr.dime.pixdim(2:4)=nii.hdr.dime.pixdim(2:4).*factor;
prefix=strtok2(fname,'.');
prefix=strtok2(prefix,'.');

save_untouch_niigz(nii,[prefix,'_ds'],1);