function Look_Lock

TR=1;
T1=2;
N=100;

M(1)=-0.5;
t(1)=0;
fa=30;
for i=2:N
  
  [tmp,ttmp]=relax(M(end),TR,T1);
  M(end+1:end+20)=tmp;
  t(end+1:end+20)=t(end)+ttmp;
  M(end+1)=M(end)*cos(fa*pi/180);
  t(end+1)=ttmp(end)+t(end);
  
end


figure;plot(t,M);
y0=ss_GEEPI(fa,TR,0,T1,100)/sin(fa*pi/180);
yaprox=y0+(M(1)-y0)*exp(-t/T1-t*log(cos(fa*pi/180))/TR);

hold on;
plot(t,yaprox,'k-');

y=y0+(M(1)-y0)*exp(-t/T1).*cos(fa*pi/180).^floor(t/TR);


plot(t,y,'g-');


function [s,t]=relax(m0,delay,T1)

t=linspace(0,delay,20);
s=1-(1-m0)*exp(-t/T1);

