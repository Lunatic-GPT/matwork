

L=1;
dk=2*pi/2/L;


N_arr=[64,96,128,256];
snr0=1000;
int=zeros(1,4);
for i=1:length(N_arr)
    N=N_arr(i);
k=(-N/2+1:N/2)*dk;
s=sin(L*k/2)./k;

s(isnan(s))=L/2;

s=sum(s);

nois=L/2/snr0*sqrt(N);

snr(i)=s/nois;
int(i)=s;
end



figure;plot(N_arr,snr.*sqrt(N_arr));