function res = ifft1c(x,dim)

% res = fft1c(x)
% 
% orthonormal forward 1D FFT
%


n=size(x,dim);
shft=zeros(1,5);
shft(dim)=ceil(n/2);

x=circshift(x,shft);

fx=ifft(x,[],dim);

fx=circshift(fx,shft);

res = sqrt(n)*fx;



