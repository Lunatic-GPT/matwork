function c=trilinear_interpolation(xyz_new,xyz,v)
% xyz_new: 3*1
% xyz: 3*2
% v: 2*2*2; (x*y*z)

xyzd=zeros(1,3);

for i=1:3
xyzd(i)=(xyz_new(i)-xyz(i,1))/(xyz(i,2)-xyz(i,1));

end


c00=v(1,1,1)*(1-xyzd(1))+v(2,1,1)*xyzd(1);
c01=v(1,1,2)*(1-xyzd(1))+v(2,1,2)*xyzd(1);
c10=v(1,2,1)*(1-xyzd(1))+v(2,2,1)*xyzd(1);
c11=v(1,2,2)*(1-xyzd(1))+v(2,2,2)*xyzd(1);

c0=c00*(1-xyzd(2))+c10*xyzd(2);
c1=c01*(1-xyzd(2))+c11*xyzd(2);

c=c0*(1-xyzd(3))+c1*xyzd(3);
