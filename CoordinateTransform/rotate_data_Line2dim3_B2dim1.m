function [mg,p,deg,cimgp,fov_center]=rotate_data_Line2dim3_B2dim1(p1,p2,B0_dir,vox_size,data,crop_size,interp,vox_size_new,prefix,use_myAffine)

% p1, p2: 1*3 position (i,j,k; not x,y,z) of the end points of a line; the line defines the z
% direction after rotate
% B0: 1*3 direction of B0; in the rotated data, B0 will be along the
% first dimension unless B0 is parallel to p1 and p2
% vox_size: 1*3
% crop_size: crop data to crop_size; before interpolate. the new FOV center
% will be the center of the line p1-p2.
% interp: interpolation factor

%%
if ~exist('use_myAffine','var')
use_myAffine=true;
end

if numel(interp)==1
   interp=interp*ones(1,3); 
end
d=(p1-p2).*vox_size;
deg=acos(dot(d,B0_dir)/sos(d,2))*180/pi;


if use_myAffine
    m=vec2Newz_B0InNewxz(d,B0_dir);
   
    m2=eye(4);
    m2(1:3,1:3)=m;
else
    m=vec2Newz_B0InNewxz(d,B0_dir);
    m=m.*repmat(vox_size(1:3)./interp,[3,1]);
    m2=eye(4);
    m2(1:3,1:3)=m;
end

c=round((p1+p2)/2);

cimg_tmp=data;

cimg=zeros([crop_size(:)',size(data,4)]);

ind=[c(:)-crop_size(:)/2,c(:)+crop_size(:)/2-1];

[ind,inds]=index_within_limit(ind,size(cimg_tmp(:,:,:,1)));


cimg(inds(1,1):inds(1,2),inds(2,1):inds(2,2),inds(3,1):inds(3,2),:)=cimg_tmp(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2),:);


fcimg=fft2c(cimg);
fcimg=fft1c(fcimg,3);

fcimg2=zpad(fcimg,crop_size(1)*interp(1),crop_size(2)*interp(2),crop_size(3)*interp(3),size(fcimg,4));

cimg2=ifft2c(fcimg2);
cimg2=ifft1c(cimg2,3);

for i=1:size(cimg2,4)
    if ~use_myAffine
     cimgp(:,:,:,i)=affine(cimg2(:,:,:,i),m2,vox_size_new.*[1,1,1],0);
      fov_center=ceil((size(cimgp(:,:,:,i))+1)/2);
    else
     [cimgp(:,:,:,i),fov_center]=myAffine(cimg2(:,:,:,i),m2,vox_size./interp,vox_size_new.*[1,1,1]); %added interp on 5/18/2018
    end
end

mg=abs(cimgp);
p=angle(cimgp);

if nargout==0
    pos_str1=sprintf('%d_',p1);
    pos_str2=sprintf('%d_',p2);
    
    if ~exist('prefix','var')
        save(['affine_',pos_str1(1:end-1),'-',pos_str2(1:end-1)],'mg', 'p','deg');
    else
        save(prefix,'mg', 'p','deg');  
    end
end


%%










