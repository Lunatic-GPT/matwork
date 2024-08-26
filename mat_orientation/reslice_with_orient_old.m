function v=reslice_with_orient_old(sz_ref,o_ref,o_in,d_in,use_nn)
% 
% d1: the underlay
% d2: the overlay; 
% sample the overlay data to match underlay grid
% same result as v=reslice_with_orient but less efficient;

if ~exist('use_nn','var')
    use_nn=true;
end

d1=o_ref;
d2=o_in;

if size(d_in,3)==1
    d_in=repmat(d_in,[1,1,2]);
    d2.pos=center2pos(d2.voxsize,d2.rotmat,size(d_in),d2.center);       
end

[xx2,yy2,zz2]=get_grid(size(d_in),d2.voxsize,eye(3),d2.rotmat'*d2.pos);

[xx,yy,zz]=get_grid(sz_ref,d1.voxsize,d2.rotmat'*d1.rotmat,d2.rotmat'*d1.pos);

method='linear';
if use_nn
    method='nearest';
end

xx3=permute(xx2,[2,1,3]);
yy3=permute(yy2,[2,1,3]);
zz3=permute(zz2,[2,1,3]);
d_in2=permute(d_in,[2,1,3,4]);

for i=1:size(d_in,4)

 v(:,:,:,i) = interp3(xx3,yy3,zz3,d_in2(:,:,:,i),xx,yy,zz,method);


end
fprintf('underlay z: %3.2f - %3.2f; overlay z: %3.2f - %3.2f\n',min_max(zz2(:)),min_max(zz(:)));    

     
      
function [xx,yy,zz]=get_grid(sz,voxsize,rotmat,pos)
        
if length(sz)==2
    sz(3)=1;
end


%[x,y,z]=meshgrid((0:sz(2)-1)*voxsize(2),(0:sz(1)-1)*voxsize(1),(0:sz(3)-1)*voxsize(3));


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

