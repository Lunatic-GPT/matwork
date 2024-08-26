function res=transform_matrix_rotation(theta,phi)
% theta goes from positive z to positive x
% phi goes from positive x to positive y
% right-handed coordinate system
% angles in degrees
% rotate coordinates first by theta around y; then by phi around z;
% matrix should be on the left, i.e. vec_new=A*vec;

theta=theta*pi/180;
phi=phi*pi/180;

m1=eye(3);

m1([1,3],[1,3])=[cos(theta),sin(theta);-sin(theta),cos(theta)];

m2=eye(3);

m2(1:2,1:2)=[cos(phi),-sin(phi);sin(phi),cos(phi)];

res=m2*m1;