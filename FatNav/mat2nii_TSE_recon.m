function mat2nii_TSE_recon(fmat,fpro,prefix)

pos=zeros(1,3);
fov=zeros(1,3);
norm=zeros(1,3);

if isa(fmat,'char')
  load(fmat); %should contain res
end


pos(1) = readsPar(fpro,'asSlice[0].sPosition.dSag');
pos(2)= readsPar(fpro,'asSlice[0].sPosition.dCor');
pos(3)=readsPar(fpro,'asSlice[0].sPosition.dTra');
fov(1) = readsPar(fpro,'asSlice[0].dReadoutFOV');
fov(2) = readsPar(fpro,'asSlice[0].dPhaseFOV');
fov(3) = readsPar(fpro,'asSlice[0].dThickness');
nslab=readsPar(fpro,'sKSpace.lImagesPerSlab'); 

sz=size(res);
sz(3)=nslab;
o.voxsize=fov./sz;
o.center=pos;

norm(1)=readsPar(fpro,'asSlice[0].sNormal.dSag');
norm(2)=readsPar(fpro,'asSlice[0].sNormal.dCor');
norm(3)=readsPar(fpro,'asSlice[0].sNormal.dTra');
dphi=readsPar(fpro,'asSlice[0].dInPlaneRot');

rotmat=get_rotmat_old(norm,dphi); %old sequence

o.rotmat=rotmat;


o.orient='SPR';
o.sz=size(res);
o=center2pos_o(o,o.center);

nii=make_nii(res); 
nii=nii_from_orient(nii,o);
save_untouch_niigz(nii,prefix);  



function rotmat=get_rotmat_old(norm,dphi)

rotmat=NormInplaneRot2Rotmat(norm,dphi);
rotmat(:,[1,2,3])=-rotmat(:,[1,2,3]); %Dichter data




