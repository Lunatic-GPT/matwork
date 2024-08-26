fov=128;

k = -(-64:63)*2*pi/fov;

dx = 1;

a=zeros(1,128);
a(32:96)=1;

fa=fft1c(a,2);

fa2=fa.*exp(-1i*k*dx);

a2=ifft1c(fa2,2);

figure;plot(a);
hold on;plot(abs(a2),'r');

