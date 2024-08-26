clear all;
N=16;   %number of phase encoding steps.
n=-N/2+1:N/2;
x=linspace(0,1,N*10);

for i=1:length(x)
    
k(i)=sum(exp(1i*x(i)*n*2*pi));
end

k=fftshift(k);
figure;plot((x-0.5)*N,real(k),'k-');