function applyXfm_flirt(in,omat,ref,use_nn,prefix)
% apply deformation field from flirt as saved in omat;
% applyDeformationField(in,omat,ref)
% in: input volume
% ref: reference volume
% omat: transformation matrix;


%pvsmask='demon/T2_PVS34_flirt.nii.gz';

%pvsmask='mask_koji_ver2/Prob_Itr1_mask_PVS01_koji2.nii';

%omat='demon/omat_PVS34.1D';

if ~exist('use_nn','var')
    use_nn=false;
end
if isa(omat,'char')
omat=load(omat);
end

if isa(in,'char')
 nii=load_untouch_niigz(in);
else
    nii=in;
end

nii_ref=load_untouch_niigz(ref);

xyz=get_coord_nii_method1(nii);
xyz_ref=get_coord_nii_method1(nii_ref);
xyz2=inverse_transform_coord(xyz_ref,omat(1:3,:));


d=fix_orientation(nii.img,nii,nii_ref);

img_new=do_interp(xyz,d,xyz2,use_nn);  %flip because data orientation different between 

nii2=nii_ref;
nii2.img=img_new;
nii2.hdr.dime.dim(6)=1;

if ~exist('prefix','var')
  prefix=strtok(in,'.');
  prefix=[prefix,'_dfm'];
end

save_untouch_niigz(nii2,prefix);

function d=fix_orientation(d,nii,nii_ref)

per=zeros(1,3);
fl=zeros(1,3);

m=[nii.hdr.hist.srow_x;nii.hdr.hist.srow_y;nii.hdr.hist.srow_z];
mref=[nii_ref.hdr.hist.srow_x;nii_ref.hdr.hist.srow_y;nii_ref.hdr.hist.srow_z];

for i=1:3  %reference
    for j=1:3
     
        costh(j)=sum(mref(i,1:3).*m(j,1:3))/sos(mref(i,1:3))/sos(m(j,1:3));
     
    end
    
        [tmp,per(i)]=max(abs(costh));
        
        fl(i)=sign(costh(per(i)));
        
end

d=permute(d,per);
for i=1:3
    if fl(i)<0
      d=flip(d,i);
    end
end

function xyz=transform_coord(xyz,m)

m2=reshape(m(:,1:3)',[1,1,1,3,3]);
pos=reshape(m(:,4),[1,1,1,3]);
xyz=squeeze(sum(m2.*xyz,4))+pos;


function xyz=inverse_transform_coord(xyz,m)

im=inv(m(:,1:3));

im2=reshape(im(:,1:3)',[1,1,1,3,3]);
pos=reshape(m(:,4),[1,1,1,3]);

xyz=squeeze(sum(im2.*(xyz-pos),4));


function res=do_interp(xyz0,d,xyznew,use_nn)


method='linear';
if use_nn
    method='nearest';
end

xx2=xyz0(:,:,:,1);
yy2=xyz0(:,:,:,2);
zz2=xyz0(:,:,:,3);

xx=xyznew(:,:,:,1);
yy=xyznew(:,:,:,2);
zz=xyznew(:,:,:,3);

res = interp3(permute(xx2,[2,1,3]),permute(yy2,[2,1,3]),permute(zz2,[2,1,3]),permute(d,[2,1,3,4]),xx,yy,zz,method);


function xyz=get_coord_nii_method3(nii)
        %%

m=[nii.hdr.hist.srow_x;nii.hdr.hist.srow_y;nii.hdr.hist.srow_z];

d=nii.img;
ijk=meshgrid2(0:size(d,1)-1,0:size(d,2)-1,0:size(d,3)-1);

xyz=transform_coord(ijk,m);

function xyz=get_coord_nii_method2(nii)
        %%
     if nii.hdr.hist.qform_code==0
         error('cannot be used when qform code is zero');
     end
voxsize=nii.hdr.dime.pixdim;
b=nii.hdr.hist.quatern_b;
c=nii.hdr.hist.quatern_c;
d=nii.hdr.hist.quatern_d;
a=sqrt(1-b^2-c^2-d^2);

R=[a^2+b^2-c^2-d^2,2*(b*c-a*d),2*(b*d+a*c);...
    2*(b*c+a*d),  a^2+c^2-b^2-d^2, 2*(c*d-a*b);...
    2*(b*d-a*c), 2*(c*d+a*b), a^2+d^2-b^2-c^2];
qoffset=[nii.hdr.hist.qoffset_x,nii.hdr.hist.qoffset_x,nii.hdr.hist.qoffset_x]';

voxsize=[voxsize(2),voxsize(3),voxsize(4)*voxsize(1)];
R=[R.*voxsize,qoffset];


ijk=meshgrid2(0:size(d,1)-1,0:size(d,2)-1,0:size(d,3)-1);

xyz=transform_coord(ijk,R);



function xyz=get_coord_nii_method1(nii)
        %%
        
        d=nii.img;
voxsize=nii.hdr.dime.pixdim;

ijk=meshgrid2(0:size(d,1)-1,0:size(d,2)-1,0:size(d,3)-1);

voxsize=reshape(voxsize(2:4),[1,1,1,3]);
xyz=ijk.*voxsize;



