function orient=get_orient(fname)


if strcmp(fname(end-2:end),'.gz') 
nii=load_untouch_niigz(fname);
orient=get_orient_from_nii(nii);

elseif strcmp(fname(end-3:end),'.nii')
nii=load_untouch_niigz(fname);
orient=get_orient_from_nii(nii);

else


orient=load(fname,'voxsize','center','rotmat','pos','orient');
    
    
end

