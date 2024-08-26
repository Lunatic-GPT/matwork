trj=readTraj('20');
x=squeeze(trj(1,:,:));
y=squeeze(trj(2,:,:));
z=squeeze(trj(3,:,:));

r=sqrt(x.^2+y.^2);
r2=sqrt(r.^2+z.^2);

theta=atan2(r,z);

phi=atan2(y,x);

figure;hist(phi(end,:));
figure;hist(theta(end,:),0:0.1:3.14159);
