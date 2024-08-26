function nii=nii_permute_data(nii,dim_permute)
% dim_permute should have the same formate as the second parameter for
% permute

   prefix=[];
if isa(nii,'char')
    prefix=strtok2(nii,'.');
    nii=load_untouch_nii(nii);
end

orient=get_orient_from_nii(nii);

orient.rotmat=orient.rotmat(:,dim_permute);


%orient.pos=orient.pos(dim_permute);

orient.voxsize=orient.voxsize(dim_permute);
nii.img=permute(nii.img,dim_permute);
nii.hdr.dime.dim(2:4)=size(nii.img);
nii.hdr.dime.pixdim(2:4)=nii.hdr.dime.pixdim(dim_permute+1);

nii=nii_from_orient(nii,orient);

    
if ~isempty(prefix)
   save_untouch_nii(nii,sprintf('%s_pmt_%d%d%d',prefix,dim_permute));
    
end
