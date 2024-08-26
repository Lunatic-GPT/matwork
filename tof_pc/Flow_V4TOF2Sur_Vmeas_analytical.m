function [f,vmean]=Flow_V4TOF2Sur_Vmeas_analytical(v,tof,r,thk,TR,fa,T1)
% [f,vmean]=FlowVelocity4VmeasTOF_analytical(v,tof,r,thk,TR,fa,T1)
% T1 is the flip angle for the flowing spins
% v: velocity in cm/s
% tof: enhancement relative to nearby voxel with only static tissue
% r: signal ratio of flowing spin to static spin when there is no flow enhancement (i.e. v=0)
% thk: thickness in cm
% TR: tr in s
% fa: flip angle in degrees
% T1 in s
% f: v*pvf (unit cm/s)
% vmean: mean velocity (cm/s)
% Note: The function only works when flow enhancement is proportional to
% velocity which is only true for low flip angles.

k=slope_TOF_vs_v(thk,TR,fa,T1);

plug_flow=false;
if plug_flow
  kp=k;  
else
  kp=4/3*k;
end

a=r*(kp*r-k*r-kp);

b=r*tof-r+v*r*k*tof;

c=-tof*(tof-1)*v;

if (tof*(1-v*k)-1)<0
 f=(-b+sqrt(b^2-4*a*c))/2/a;
else
 f=(-b-sqrt(b^2-4*a*c))/2/a;
end
x=f;
a2=kp*x*r;
b2=(r*x-v-v*k*r*x);
c2=-v*(r-1)*x;

vmean=(-b2+sqrt(b2^2-4*a2*c2))/2/a2;

function k=slope_TOF_vs_v(thk,TR,fa,T1)
v=linspace(0.01,0.3,20);
v=v(:);
nv=length(v);


for i=1:length(v)
  sart(i)= TOF_signalIntensity(v(i),thk,TR,fa,T1);
end

xmat=[ones(nv,1),v(1:nv)];

b=xmat\sart(:);
% 
% 
% figure;plot(v,sart,'o-');
% hold on;
% 
% plot(v,xmat*b,'r-');
k=b(2)/b(1);




