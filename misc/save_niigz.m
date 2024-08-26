function save_niigz(nii,prefix)

save_nii(nii,[prefix,'.nii']);
gzip([prefix,'.nii']);
delete([prefix,'.nii']);