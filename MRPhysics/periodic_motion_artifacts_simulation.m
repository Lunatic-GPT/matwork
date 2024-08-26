a=zeros(1,1024);
a(256)=1;
fa=ifft1c(a,2);

dp=phase(fa(1))-phase(fa(2));
b=fa.*exp(1i*(-512:511)*dp/128.*sin((-512:511)/1024*20*2*pi));

fb=fft1c(b,2);

figure;plot(abs(fb));

