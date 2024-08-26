
t=1:121;
a=exp(-t/121*log(2));

figure;plot(a);

ae=[a*0.5,a,a*0.5];

figure;plot(ae);
%%
hold on;plot(0:362,abs(ifft1c(ae,2))/3.7,'g');

%%
gamma=42.58e6;  %s-1/T
fov=0.21;   % m

gamma=gamma*1e-9;  %us-1/mT

 mom=1/gamma/fov*420/3; %us*mT/cm
 
 
 
 

