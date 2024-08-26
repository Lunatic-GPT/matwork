function v2=reslice(v,m,vox,center,sz2,vox2,twoD,use_nn)
% v2=reslice(v,m,vox,center_new,sz2,vox2)
% where v is the data
% m is the rotation matrix; to transform from the old grid coordinate 
% to the new grid coordinate system;
% x,y,z are the 1, 2, and 3 dimensions;

% vox is the voxel size for v
% center: new center
% sz2: matrix size of resliced data
% vox2: voxel size of resliced data
% the center of the old data set is assumed to be at voxel (N+1)/2 (odd
% dim) or N/2+1 (even dim) and equal to 0;
% twoD: whether to interpolate along the third dimension
% use_nn:  whether to use nearest-neighbor interpolation; otherwise linear.

% m(1,2)=-m(1,2);
% m(2,1)=-m(2,1);
% m(3,2)=-m(3,2);
%m(2,3)=-m(2,3);

if length(sz2)==2
    sz2(3)=1;
end

if twoD
    
   if vox(3)~=vox2(3)  || sz2(3)~=size(v,3)
       error('voxel size of matrix dimension mismatch for 2D reslice');
   end
   
       
end
x=get_vec(size(v,1))*vox(1);  %pos for the original data
y=get_vec(size(v,2))*vox(2);
z=get_vec(size(v,3))*vox(3);

[xx,yy,zz]=meshgrid(x,y,z);

x2=get_vec(sz2(1))*vox2(1)+center(1);   %pos for the resampled data
y2=get_vec(sz2(2))*vox2(2)+center(2);
z2=get_vec(sz2(3))*vox2(3)+center(3);

xyz2=transform_grid(inv(m),x2,y2,z2);

method='linear';
if use_nn
    method='nearest';
end

v=permute(v,[2,1,3,4]);
if size(v,3)==1  || twoD
    v2=zeros(sz2);
    
    for j=1:size(v,4)
        for i=1:size(v,3)
            v2(:,:,i,j) = interp2(xx(:,:,i),yy(:,:,i),v(:,:,i,j),squeeze(xyz2(1,:,:,i)),squeeze(xyz2(2,:,:,i)),method);
        end
    end
    
else
    v2 = interp3(xx,yy,zz,v,squeeze(xyz2(1,:,:,:)),squeeze(xyz2(2,:,:,:)),squeeze(xyz2(3,:,:,:)),method);
end


          
function v2=get_vec(v)

if mod(v,2)==1
    v2=-(v-1)/2:(v-1)/2;
else
     v2=-v/2:v/2-1;
end

function xx=transform_grid(m,x,y,z)


xx=zeros([3,length(x),length(y),length(z)]);

for i=1:length(x)
    for j=1:length(y)
        for k=1:length(z)
          xx(:,i,j,k)=[x(i),y(j),z(k)]*m;
          
        end
    end
end

