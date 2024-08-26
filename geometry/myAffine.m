function [new_img,orig] = myAffine(old_img, M, elem_size, new_elem_size)

% M is the matrix that rotates the old image coordinates (i.e. in step of elem_size)
% to the new image frame (in step of new_elem_size)


cm=floor(size(old_img)/2)+1;

x0=(1:size(old_img,1))-cm(1);
y0=(1:size(old_img,2))-cm(2);
z0=(1:size(old_img,3))-cm(3);
x0=x0*elem_size(1);
y0=y0*elem_size(2);
z0=z0*elem_size(3);

[xx,yy,zz]=meshgrid2(x0,y0,z0);

[xx2,yy2,zz2]=coord_transform(xx,yy,zz,M);

xl=min_max(xx2(:));
yl=min_max(yy2(:));
zl=min_max(zz2(:));

xq=round(xl(1)/new_elem_size(1)):round(xl(2)/new_elem_size(1));
yq=round(yl(1)/new_elem_size(2)):round(yl(2)/new_elem_size(2));
zq=round(zl(1)/new_elem_size(3)):round(zl(2)/new_elem_size(3));

xq=xq*new_elem_size(1);
yq=yq*new_elem_size(2);
zq=zq*new_elem_size(3);


[xq,yq,zq]=meshgrid2(xq,yq,zq);

[xq,yq,zq]=coord_transform(xq,yq,zq,inv(M));

new_img=interp3(yy,xx,zz,old_img,yq,xq,zq);

r=sqrt(xq.^2+yq.^2+zq.^2);

[tmp,orig]=min(r(:));
orig=ind2subb(size(r),orig);



function [xx2,yy2,zz2]=coord_transform(xx,yy,zz,M)



xx2=shiftdim(xx,-2);
yy2=shiftdim(yy,-2);
zz2=shiftdim(zz,-2);

xyz=cat(2,xx2,yy2,zz2);
xyz(1,4,:,:,:)=1;

                  rls=version('-release');
                  rls=str2num(rls(1:4));
                  if rls<=2013
                      sz=size(xyz);
                      sz(1:2)=1;
                     M=repmat(M,sz);
                     xyz=repmat(xyz,[4,1,1,1,1]);
                     xyz_new=squeeze(sum(M.*xyz,2));
                     
                  else
                     xyz_new=squeeze(sum(M.*xyz,2));
                  end
xx2=squeeze(xyz_new(1,:,:,:));
yy2=squeeze(xyz_new(2,:,:,:));
zz2=squeeze(xyz_new(3,:,:,:));
