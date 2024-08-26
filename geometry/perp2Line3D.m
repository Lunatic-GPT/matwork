
function res=perp2Line3D(v)

% generate two mutually perpendicular vectors that are both perp to v
% res: 2*3
v=v(:)';

v=v(:)'/sqrt(sum(v.*v));

[~,i]=max(v);

t=ones(1,3);

t(i)=0;

t=t-v*t'*v;
t=t/sqrt(sum(t.*t));

t2=cross(t,v);
t2=t2/sqrt(sum(t2.*t2));

res=[t;t2];


