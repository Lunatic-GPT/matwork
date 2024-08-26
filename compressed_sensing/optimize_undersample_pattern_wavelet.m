XOP = Wavelet_ISU([448,448]);
m=vdPoisMex(448,448,1,1,2,2,[0,0],0,4,0);

x=m*0;
x(end/2,end/2)=1;

fx=fft2c(x);
fx2=fx.*m;
x2=ifft2c(fx2);

wx2=XOP*x2;

wx=XOP*x;

figure;imshow(abs(wx2),[0,0.002]);
