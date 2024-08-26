function dcm2nii_mosaic_manual(fname,pro,main_inplane)
% tested for axial with inplane = 0 for both main and nav;

if isa(main_inplane,'char')
  main_inplane=str2num(main_inplane);
end

in=dicominfo(fname);
a=in.ImageOrientationPatient;
b=cross(a(1:3),a(4:6));
o.rotmat=[a(1:3),a(4:6),-b];


o.voxsize=in.PixelSpacing(1)*[1,1,1];

x=readsPar(pro,'sNavigatorArray.asElm[0].sCuboid.sPosition.dSag');
y=readsPar(pro,'sNavigatorArray.asElm[0].sCuboid.sPosition.dCor');
z=readsPar(pro,'sNavigatorArray.asElm[0].sCuboid.sPosition.dTra');


d=dicomread(fname);


d=mosaic_split(d,8);
d=flip(d,3);
d=flip(d,2);
d=flip(d,1);

if main_inplane==0

center=[y,x,-z];

elseif main_inplane==90
    center=[0,x,y]; %still not correct;
   % d=flip(d,1);
end


nii=make_nii(d);

o.pos=center2pos(o.voxsize,o.rotmat,size(d),center);
nii=nii_from_orient(nii,o);


prefix=strtok(filename(fname));
save_untouch_nii(nii,[prefix,'.nii']);

