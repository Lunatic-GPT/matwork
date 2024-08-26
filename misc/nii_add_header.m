function nii_add_header(orig,parent)

try 
nii=load_nii(parent);
catch
nii=load_untouch_nii(parent);
end

nii_orig=load_nii(orig);
nii.img=nii_orig.img;

nii.hdr.dime.datatype=nii_orig.hdr.dime.datatype;
nii.hdr.dime.bitpix=nii_orig.hdr.dime.bitpix;
    
[prefix,suf]=strtok(orig,'.');

try
save_nii(nii,[prefix,'2',suf]);
catch
save_untouch_nii(nii,[prefix,'2',suf]);
    
end