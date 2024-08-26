function applyDeformationField(in,omat,dfield,boxsize,mask,prefix,debug)
% obsolete; use applyDeformationField_general instead
%T2 to atlas: 
% Not sure why need a different (from applyDeformationField_atlas2T2) one if transform from T2 to Atlas;

% demon appears to adopt the world coordinate system;

if ~exist('boxsize','var')
    boxsize=6;
end

if ~exist('debug','var')
    debug=0;
end

if ~exist('prefix','var')
prefix=[strtok(filename(in),'.'),'_Eve'];
end


if isa(omat,'char')
omat=load(omat);
omat=omat(1:3,:);
end
nii=load_untouch_niigz(in);
nii.img=nii.img>0;
nii_df=load_untouch_niigz(dfield);

%xyz=get_coord_nii_method1(nii);
%xyz2=transform_coord(xyz,omat(1:3,:)); %do it inside transform demon;
%d=nii.img;
d=fix_orientation(nii.img,nii,nii_df);
 
voxsize=nii_df.hdr.dime.pixdim(2:4);

voxsize_d=nii.hdr.dime.pixdim(2:4);


for i=1:3
  df(:,:,:,i)=nii_df.img(:,:,:,i)*voxsize(i);
end

if ~exist('mask','var') || isempty(mask)
    mask=ones(size(df(:,:,:,1)));
else
    mask=ri(mask);
end
  mask=ones(size(df(:,:,:,1)));
[df,mask]=fix_orientation_coord(df,mask>0,nii_df);

%img_new=transform_demon_interp(xyz,d,df,omat,voxsize);

img_new=transform_demon_average(d,voxsize_d,df,omat,voxsize,boxsize,mask,debug);
%img_new=do_average(xyz3,d,voxsize_df,size(nii_df.img(:,:,:,1)));  % calculate average over voxels that resides in the voxels defined by xyz_ref;   

%img_new=count_overlap(xyz3,d,xyz_ref);  % count voxels that resides in the voxels defined by xyz_ref;   
%img_new=fix_orientation(img_new,nii,nii_df);

nii2=nii_df;
nii2.img=img_new;
nii2.hdr.dime.dim(6)=1;


save_untouch_nii(nii2,prefix);

function res=transform_demon_average(d,voxsize_d,df,omat,voxsize,boxsize,mask,debug)
% deformatoin field in units of voxsize;

dfx=df(:,:,:,1);
dfy=df(:,:,:,2);
dfz=df(:,:,:,3);

sz=size(d);
res=0*df(:,:,:,1);

for i=1:size(df,1)
        if debug && i~=round(size(df,1)/2)
           continue;
        end
    
    disp(i);
    for j=1:size(df,2)
        for k=1:size(df,3)
          %  tic;
          if mask(i,j,k)==0
              continue;
          end
            xyz0=[i-1,j-1,k-1].*voxsize(:)';
            xyz_tmp=xyz0'+[dfx(i,j,k),dfy(i,j,k),dfz(i,j,k)]';
            xyz_new=inverse_transform_coord2D(xyz_tmp,omat);
            
            ijk_d=round(xyz_new'./voxsize_d)+1;
            
           [i1,j1,k1]=get_ijk(ijk_d,boxsize,sz);
           
           if ~isempty(i1)&&~isempty(j1)&&~isempty(k1)
            res(i,j,k)=mean(vec(d(i1,j1,k1)));
           end
            
          %  toc;
        end
    end
end

function [i1,j1,k1]=get_ijk(ijk_d,boxsize,sz)

 i1=ijk_d(1)-boxsize/2:ijk_d(1)+boxsize/2;
            j1=ijk_d(2)-boxsize/2:ijk_d(2)+boxsize/2;
            k1=ijk_d(3)-boxsize/2:ijk_d(3)+boxsize/2;
            
            i1(i1<1)=[];
            i1(i1>sz(1))=[];
            
            j1(j1<1)=[];
            j1(j1>sz(2))=[];
           
            k1(k1<1)=[];
            k1(k1>sz(3))=[];

 
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


function [df,mask]=fix_orientation_coord(df,mask,nii)

per=zeros(1,3);
fl=zeros(1,3);

m=[nii.hdr.hist.srow_x;nii.hdr.hist.srow_y;nii.hdr.hist.srow_z];
mref=[-1,0,0,0;0,-1,0,0;0,0,1,0];  % Demon uses dicom convention; 

for i=1:3  %reference
    for j=1:3
     
        costh(j)=sum(mref(i,1:3).*m(j,1:3))/sos(mref(i,1:3))/sos(m(j,1:3));
     
    end
    
        [tmp,per(i)]=max(abs(costh));
        
        fl(i)=sign(costh(per(i)));
        
end

df=permute(df,[per,4,5]);
for i=1:3
    if fl(i)<0
      df(:,:,:,i)=-df(:,:,:,i);
    end
end

mask=permute(mask,[per,4,5]);

function xyz=transform_coord(xyz,m)

m2=reshape(m(:,1:3)',[1,1,1,3,3]);
pos=reshape(m(:,4),[1,1,1,3]);
xyz=squeeze(sum(m2.*xyz,4))+pos;


function xyz=inverse_transform_coord(xyz,m)


im=inv(m(:,1:3));

im2=reshape(im(:,1:3)',[1,1,1,3,3]);
pos=reshape(m(:,4),[1,1,1,3]);

xyz=squeeze(sum(im2.*(xyz-pos),4));


function xyz=inverse_transform_coord2D(xyz,m)


pos=m(:,4);

xyz=m(:,1:3)\(xyz-pos);

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


function res=transform_demon_interp(xyz,d,df,omat,voxsize,interp)
% deformatoin field in units of voxsize;
voxsize=squeeze(voxsize)';
xyz_new=df*0;
dfx=df(:,:,:,1);
dfy=df(:,:,:,2);
dfz=df(:,:,:,3);

for i=1:size(df,1)
    disp(i);
    for j=1:size(df,2)
        for k=1:size(df,3)
          %  tic;
            xyz0=[i-1,j-1,k-1].*voxsize;
            xyz_tmp=xyz0'+[dfx(i,j,k),dfy(i,j,k),dfz(i,j,k)]';
            xyz_new(i,j,k,:)=inverse_transform_coord2D(xyz_tmp,omat);
          %  toc;
        end
    end
end

 res=do_interp(xyz,d,xyz_new,true);


