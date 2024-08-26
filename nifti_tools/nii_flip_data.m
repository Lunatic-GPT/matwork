function nii=nii_flip_data(nii,dim_flip)
% flip data but keep their physical coordinates the same
% 


prefix=[];
if isa(nii,'char')
    prefix=strtok2(nii,'.');
    nii=load_untouch_nii(nii);
end

orient=get_orient_from_nii(nii);

ind=[0,0,0];

sz=size(nii.img);

rotmat0=orient.rotmat;
orient.rotmat(:,dim_flip)=-orient.rotmat(:,dim_flip);

nii.img=flip(nii.img,dim_flip);
ind(dim_flip)=sz(dim_flip)-1;

orient.pos = orient.pos+(rotmat0.*orient.voxsize)*ind';
nii=nii_from_orient(nii,orient);

if ~isempty(prefix)
   save_untouch_nii(nii,sprintf('%s_flip%d',prefix,dim_flip));  
end
    
    

