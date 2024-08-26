function res=dist2plane(v,center,norm)
% v=n*3;
% center: 1*3: a point on the plane
% norm=1*3: the normal direction of the plane
% has a sign.

res=sum((v-center(:)').*norm(:)',2);