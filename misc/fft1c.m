function res = fft1c(x,dim)

% res = fft1c(x,dim)
% 
% orthonormal forward 1D FFT
%
% (c) Michael Lustig 2005

n=size(x,dim);
shft=zeros(1,4);
shft(dim)=-ceil(n/2);

x=circshift(x,shft);

fx=fft(x,[],dim);

fx=circshift(fx,shft);

res = 1/sqrt(n)*fx;



