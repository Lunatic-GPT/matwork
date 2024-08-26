function var2niigz(d,parent_nii,prefix)


if ~exist('parent_nii','var')
   save_nii(make_nii(d),sprintf('%s.nii',prefix));
else
    
    if strcmp(parent_nii(end-2:end),'.gz')
     nii=load_untouch_niigz(parent_nii);
    else
     nii=load_untouch_nii(parent_nii);
    end
  %  nii.img=flip(flip(permute(d,[2,1,3,4]),2),1);
   
    nii.img=d;
    save_untouch_niigz(nii,prefix);  % the orientation of the produced nii file differs from the original one. not good!  Ok if use as underlay and overlay together.
   
end
