function res=transform_matrix_rotation_arb_axis(theta,phi,deg)
% rotating by deg degrees around an axis defined by the angles theta and phi

% the following is ok too
% zp=thetaPhi2unitVec(theta,phi);
% 
% [tmp,ind]=max(zp);
% yp=ones(1,3);
% yp(ind)=0;
% 
% yp=yp/sqrt(2);
% yp=cross(yp,zp);
% yp=yp/sos(yp,2);
% xp=cross(yp,zp);
% 
% m1=[xp',yp',zp'];

m1=transform_matrix_rotation(theta,phi);

m2=transform_matrix_rotation(0,deg);

res=m1*m2/m1;

