function ice2nii_mosaic(dname,ColVec,RowVec,NormalVec,voxSize,center)

%add later:
%name=name4pat(fullfile(dname,'*.IceHead'));
%row=readxPar(name{1},'ColVec');
% pos in IceHead is the center of slice;
%% this to be read from IceHead
center=[0,-24.289,-2.6988]; %pos in IceHead gave the center of the most inferior slice.
center=[0,0,0];
ColVec=[0,1,0];
RowVec=[-1,0,0];
NormalVec=[0,0,1];
voxSize=[4,4,4];
reps=8;
%%
a=ri_ice(dname);
a=mosaic_split(a,reps);
a=permute(a,[2,1,3,4]);
a=flip(a,3);
a=flip(a,1);

%a=flip(a,2);
nii=make_nii(a);

orient.voxsize=voxSize;
 orient.rotmat=[ColVec(:),RowVec(:),NormalVec(:)];
 orient.pos=center2pos(voxSize,orient.rotmat,size(a),center);
 
nii=nii_from_orient(nii,orient);


mat2niigz(a,'',dname,0,nii);
    
