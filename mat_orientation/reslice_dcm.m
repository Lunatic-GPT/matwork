function v=reslice_dcm(d1,d2,prefix,use_nn,nslice,d1var,d2var,d1center_new,rotmat_new)
%  v=reslice_dcm(d1,d2,prefix,use_nn,nslice,d1var,d2var,d1centerr_new,rotmat_new)
% d1: the underlay; can be .mat, dicom dir, .nii, or .niigz.
% d2: the overlay; 
% d1center_new: new center position for the d1
% rotmat_new: new rotation matrix for d1; can be obtained from
% NormInplaneRot2Rotmat;
% sample the overlay data to match underlay grid


if ~exist('use_nn','var')
   use_nn=true;
end
[d1,tmp]=strtok(d1,'$');
if ~isempty(tmp)
    d1var=tmp(2:end);
end

[d2,tmp]=strtok(d2,'$');
if ~isempty(tmp)
    d2var=tmp(2:end);
end

if ~exist('d1var','var') || isempty(d1var)
   d1var='d';
end

if ~exist('d2var','var') || isempty(d2var)
   d2var='d';
end

nii=[];
if (length(d1)>=4 && strcmp(d1(end-3:end),'.nii')) || (length(d1)>=7&&strcmp(d1(end-6:end),'.nii.gz')) 
    nii=load_untouch_niigz(d1);
end

 d1=get_data(d1);
 if exist('d1center_new','var')
   d1.center=d1center_new;
 end
 if exist('rotmat_new','var')
   d1.rotmat=rotmat_new;
 end
 
d1.pos=center2pos(d1.voxsize,d1.rotmat,size(d1.d),d1.center);
d1d=getfield(d1,d1var);
d2=get_data(d2);
d2d=getfield(d2,d2var);

d2d=double(d2d);
[xx2,yy2,zz2]=get_grid(size(d2d),d2.voxsize,eye(3),d2.rotmat'*d2.pos);

d2.rotmat=real(d2.rotmat);
if ~exist('nslice','var') || isempty(nslice)
   dnew=d1d;
   pos=d1.pos;
else
   sz=size(d1d);
   sz(3)=nslice;
   dnew=zeros(sz);
   pos=d1.pos-d1.rotmat(:,3)*ceil((nslice-size(d1d,3))/2)*d1.voxsize(3);
end
sz_dnew=size(dnew);
if length(sz_dnew)==2
    sz_dnew(3)=1;
end
[xx,yy,zz]=get_grid(sz_dnew,d1.voxsize,d2.rotmat'*d1.rotmat,d2.rotmat'*pos);
xx3=squeeze(xx2(:,1,1));
yy3=squeeze(yy2(1,:,1));
zz3=squeeze(zz2(1,1,:));

%m_nn = find_nn(xx3,yy3,zz3,xx,yy,zz);
%save m_nn m_nn

method='linear';
if use_nn
    method='nearest';
end

if size(xx2,3)>1
for i=1:size(d2d,4)
   v(:,:,:,i) = interp3(permute(xx2,[2,1,3]),permute(yy2,[2,1,3]),permute(zz2,[2,1,3]),permute(d2d(:,:,:,i),[2,1,3]),xx,yy,zz,method);
end

else
for i=1:size(d2d,4)
   v(:,:,:,i) = interp2(permute(xx2,[2,1,3]),permute(yy2,[2,1,3]),permute(d2d(:,:,:,i),[2,1,3]),xx,yy,method);
end
    
end

if exist('prefix','var') && ~isempty(prefix)
    if ~isempty(nii)
    if size(v,4)==1 && size(v,3)==1%avoid slice position error in ITKsnap when there is a single slice and 
   %  v=repmat(v,[1,1,1,2]);
    end

    nii.img=v;
       % nii=nii_from_orient(nii,d1);


        save_untouch_niigz(nii,prefix); 
    else
        
        d1=rmfield(d1,d1var);
        d1.d=v;
        save(prefix,'-struct','d1');
    end
end

function m=find_nn(xx2,yy2,zz2,xx,yy,zz)
%generate a mask for nearest-neighbors
m=zeros(length(xx2),length(yy2),length(zz2));


for i=1:size(xx,1)
       tic;
    for j=1:size(xx,2)
        for k=1:size(xx,3)
         
         
            [~,i1]=min(abs(xx(i,j,k)-xx2));
            [~,i2]=min(abs(yy(i,j,k)-yy2));
            [~,i3]=min(abs(zz(i,j,k)-zz2));
            
          m(i1,i2,i3)=1;
         
          
        end
    end
     toc;
end


function d1=get_data(d1)

if isa(d1,'char') && exist(d1,'dir')
    if ~exist([d1,'.mat'],'file')
     ri(d1,[]);  %dicom folder
    end
    d1=load(d1);
elseif strcmp(d1(end-3:end),'.mat')
    d1=load(d1);
elseif   strcmp(d1(end-2:end),'.gz')
   nii=load_untouch_niigz(d1);
   d1=get_orient_from_nii(nii);
   d1.d=nii.img;
elseif strcmp(d1(end-3:end),'.nii')
   nii=load_untouch_nii(d1); 
   d1=get_orient_from_nii(nii);
   d1.d=nii.img;
end  


function [xx,yy,zz]=get_grid(sz,voxsize,rotmat,pos)
        
if length(sz)==2
    sz(3)=1;
end
xyz=meshgrid2((0:sz(1)-1)*voxsize(1),(0:sz(2)-1)*voxsize(2),(0:sz(3)-1)*voxsize(3));

%%

rotmat=reshape(rotmat',[1,1,1,3,3]);
xyz=sum(rotmat.*xyz,4);

%%


% x=rotmat**vox(1);  %pos for the original data
% y=get_vec(size(v,2))*vox(2);
% z=get_vec(size(v,3))*vox(3);
xx=xyz(:,:,:,1,1)+pos(1);
yy=xyz(:,:,:,1,2)+pos(2);
zz=xyz(:,:,:,1,3)+pos(3);

