function applyDeformationField_general(in,omat,dfield,boxsize,prefix,omat_only,debug)
% applyDeformationField_general(in,omat,dfield,boxsize,prefix,omat_only,debug)
% in, dfield in dicom format
% omat: output from flirt; a text file containing a 4*4 matrix
% boxsize: average over how many neighboring voxels along each dimension
% dfield is the deformation field output from demon; 
% if omat_only is true, then only apply the omat.
% the code tested for aligning Eve template to PVS_R21 dataset
% demon appears to adopt the world (dicom) coordinate system;
% NOTE: make sure the data orientation (quaternion and s_row) for "in" and the image that was used
% to generate dfield and omat match.  

if ~exist('debug','var')
    debug=0;
end

if ~exist('boxsize','var')
    boxsize=1;
end

if ~exist('omat_only','var')
    omat_only=false;
end

if ~exist('prefix','var')
prefix=[strtok(filename(in),'.'),'_deform'];
end

if isa(omat,'char')
omat=load(omat);
omat=omat(1:3,:);
end
nii=load_untouch_niigz(in);
nii_df=load_untouch_niigz(dfield);

d=nii.img;
%d=fix_orientation(nii.img,nii,nii_df);
 
voxsize=nii_df.hdr.dime.pixdim(2:4);

voxsize_d=nii.hdr.dime.pixdim(2:4);

if omat_only
    nii_df.img=0*nii_df.img;
end

df=nii_df.img;
%{
fprintf('add temporary for debug - not correct\n');
for i=1:3
 df(:,:,:,1,i)=df(:,:,:,1,i)*voxsize(i);    
end
%}

%mask=ri(mask);
mask=ones(size(df(:,:,:,1)));
df=fix_orientation_coord(df,nii_df);

%

 
if do_flip_necessary(nii_df) %flirt will flip the first axis if determinant of the S matrix >0
  df=flip(df,1);
  df(:,:,:,1)=-df(:,:,:,1);
end

if do_flip_necessary(nii)
  d=flip(d,1);
end

%img_new=transform_demon_interp(xyz,d,df,omat,voxsize);

img_new=transform_demon_average(d,voxsize_d,df,omat,voxsize,boxsize,mask,debug);

if do_flip_necessary(nii_df) %flirt will flip the first axis if determinant of the S matrix >0
  img_new=flip(img_new,1); %flip back;
end

nii2=nii_df;
nii2.img=img_new;
nii2.hdr.dime.dim(6)=1;
save_untouch_nii(nii2,prefix);

function res=do_flip_necessary(nii)

orient=get_orient_from_nii(nii);
res= det(orient.rotmat(:,1:3))>0;
    



function res=transform_demon_average(d,voxsize_d,df,omat,voxsize,boxsize,mask,debug)
% deformatoin field in units of voxsize;

dfx=df(:,:,:,1);
dfy=df(:,:,:,2);
dfz=df(:,:,:,3);

sz=size(d);
res=0*df(:,:,:,1);

tic;
for i=1:size(df,1)    
    for j=1:size(df,2)
        for k= 1:size(df,3)
            
            if debug>0 && k~=debug
                continue;
            end
            
            %  tic;
            if mask(i,j,k)==0
                continue;
            end
            % flirt calculates the omat based on coordinates calculated as
            % the following line; futhermore, if the axes do not follow
            % right-hand rule, data will be flipped in the  x-axis to make
            % it so.
            % omat is defined as x_ref = omat*x_in;
            
            xyz0=[i-1,j-1,k-1].*voxsize(:)'; % coordinators of the reference
            xyz_tmp=xyz0'+[dfx(i,j,k),dfy(i,j,k),dfz(i,j,k)]'; %coordinators of the input image
            xyz_new=inverse_transform_coord2D(xyz_tmp,omat); %
            
            ijk_d=round(xyz_new'./voxsize_d)+1;
            
            
            
            [i1,j1,k1]=get_ijk(ijk_d,boxsize,sz);
            
            if ~isempty(i1)&&~isempty(j1)&&~isempty(k1)
                res(i,j,k)=mean(vec(d(i1,j1,k1)));
                
            end
            
        end
    end
    
  fprintf('Time: passed/total = (%f/%f)\n',toc, toc*size(df,1)/i);
end


function [i1,j1,k1]=get_ijk(ijk_d,boxsize,sz)

            l=ceil((boxsize-1)/2);
            i1=ijk_d(1)-l:ijk_d(1)+l;
 
            j1=ijk_d(2)-l:ijk_d(2)+l;
            k1=ijk_d(3)-l:ijk_d(3)+l;
            
            i1(i1<1)=[];
            i1(i1>sz(1))=[];
            
            j1(j1<1)=[];
            j1(j1>sz(2))=[];
           
            k1(k1<1)=[];
            k1(k1>sz(3))=[];

 


function df2=fix_orientation_coord(df,nii)


dicom_axis=zeros(1,3);  % the dicom axis ([x,y,z]) for the 1st, 2nd, and 3rd data axis.
fl=zeros(1,3);

m=[nii.hdr.hist.srow_x;nii.hdr.hist.srow_y;nii.hdr.hist.srow_z];

mref=[-1,0,0,0;0,-1,0,0;0,0,1,0]; %dicom axis;
for j=1:3  %first, second, and third data axis;
    
    
    for i=1:3
        costh(i)=sum(mref(:,i).*m(:,j))/sos(mref(:,i),1)/sos(m(:,j),1);
    end
    
        [~,dicom_axis(j)]=max(abs(costh));  
        fl(j)=sign(costh(dicom_axis(j)));
        
end
df2=df*0; % In df2, [1,2,3] corresponds to first, second, and third data axes; and pointing (i.e. positive) toward inreasing index value

for j=1:3         
      df2(:,:,:,j)=fl(j)*df(:,:,:,dicom_axis(j));
end

%mask=permute(mask,[per,4,5]);


function xyz=inverse_transform_coord2D(xyz,m)


pos=m(:,4);

xyz=m(:,1:3)\(xyz-pos);
