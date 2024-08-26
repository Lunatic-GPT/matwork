function d2=rotate_images(d,ax,phi,use_nn)
% d2=rotate_images(d,ax,phi)
% ax: 1,2, or 3
% ph: angle in degrees

an=phi*pi/180;

c=cos(an);
s=sin(an);

ind=1:3;

ind(ax)=[];

mat=eye(3);
mat(ind(1),ind(2))=s;
mat(ind(2),ind(1))=-s;
mat(ind(1),ind(1))=c;
mat(ind(2),ind(2))=c;

voxsize=1;



[xx,yy,zz]=get_grid(size(d),voxsize,eye(3));

[xxi,yyi,zzi]=get_grid(size(d),voxsize,mat);

method='linear';
if use_nn
    method='nearest';
end

d2 = interp3(permute(xx,[2,1,3]),permute(yy,[2,1,3]),permute(zz,[2,1,3]),permute(double(d),[2,1,3]),xxi,yyi,zzi,method);

d2(isnan(d2))=0;


function [xx,yy,zz]=get_grid(sz,voxsize,rotmat)

if length(sz)==2
    sz(3)=1;
end
xyz=meshgrid2(-sz(1)/2:sz(1)/2-1,-sz(2)/2:sz(2)/2-1,-sz(3)/2:sz(3)/2-1)*voxsize;

%%

rotmat=reshape(rotmat',[1,1,1,3,3]);
xyz=sum(rotmat.*xyz,4);


xx=xyz(:,:,:,1,1);
yy=xyz(:,:,:,1,2);
zz=xyz(:,:,:,1,3);

