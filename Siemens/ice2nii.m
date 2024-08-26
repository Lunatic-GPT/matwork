function ice2nii(dname,ColVec,RowVec,NormalVec,voxSize,center)

%add later:
%name=name4pat(fullfile(dname,'*.IceHead'));
%row=readxPar(name{1},'ColVec');
% pos in IceHead is the center of slice;
%% this to be read from IceHead
center=[0,-22.94,33.735];
ColVec=[0,1,0];
RowVec=[1,0,0];
NormalVec=[0,0,1];
voxSize=[1.64,1.64,1.6];

%%
a=ri_ice(dname);
nii=make_nii(a);

orient.voxsize=voxSize;
 orient.rotmat=[ColVec(:),RowVec(:),NormalVec(:)];
 orient.pos=center2pos(voxSize,orient.rotmat,size(a),center);
 
nii=nii_from_orient(nii,orient);


mat2niigz(a,'',dname,0,nii);
    
