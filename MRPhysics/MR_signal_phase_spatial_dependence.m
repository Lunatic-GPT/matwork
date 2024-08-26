function MR_signal_phase_spatial_dependence()
s=[0,0,0];

lam=3;

k=2*pi/lam;
w=10;
r1=[-1,1,0];
m=[1,0,0];
t=2;

H1=Hfield(m,k,r1)*exp(1i*w*t);

m=[conj(H1(1)),0,0];
H1x=Hfield(m,k,-r1);
disp(angle(H1x(1)));
m=[0,H1(1)*1i,0];
H1y=Hfield(m,k,-r1);

H1=H1x+H1y;


function res=Hfield(m,k,r1)

n=r1/sqrt(sum(r1.^2));
r=sqrt(sum(r1.^2));


t1=cross(cross(n,m),n)*k^2*exp(1i*k*r)/r;
t2=(3*n*dot(n,m)-m)*(1/r^3-1i*k/r^2)*exp(1i*k*r);
res=1/4/pi*(t1+t2);

