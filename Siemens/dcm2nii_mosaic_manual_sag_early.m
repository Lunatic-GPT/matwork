function dcm2nii_mosaic_manual_sag_early(fname,pro,main_inplane)
% sag; rot = 90; tested not yet passed.

if isa(main_inplane,'c')
  main_inplane=str2num(main_inplane);
end

in=dicominfo(fname);
a=in.ImageOrientationPatient;
a=round(a);
b=cross(a(1:3),a(4:6));
o.rotmat=[a(4:6),a(1:3),b];


o.voxsize=in.PixelSpacing(1)*[1,1,1];

x=readsPar(pro,'sNavigatorArray.asElm[0].sCuboid.sPosition.dSag');
y=readsPar(pro,'sNavigatorArray.asElm[0].sCuboid.sPosition.dCor');
z=readsPar(pro,'sNavigatorArray.asElm[0].sCuboid.sPosition.dTra');


d=dicomread(fname);


d=mosaic_split(d,8);
%d=permute(d,[2,1,3]);
d=flip(d,1);
%d=flip(d,2);
%d=flip(d,1);

if main_inplane==90

center=[x,z,y];

elseif main_inplane==90
 %   center=[0,x,y]; %still not correct;
   % d=flip(d,1);
end


nii=make_nii(d);

o.pos=center2pos(o.voxsize,o.rotmat,size(d),center);
nii=nii_from_orient(nii,o);


prefix=strtok(filename(fname));
save_untouch_nii(nii,[prefix,'.nii']);

