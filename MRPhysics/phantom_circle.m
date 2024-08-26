function [kdata,img] = phantom_circle( N,n_r )
%  [kdata,img] = phantom_circle( N,n_r )

img=zeros(N,N);
[x,y]=meshgrid(-N/2:N/2-1,-N/2:N/2-1);
img(x.^2+y.^2<n_r^2)=1;

kdata=fft2c(img);

