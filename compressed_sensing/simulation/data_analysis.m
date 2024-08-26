load synthesize_kdata_mask2_nih_nois0.001_bold0.010.mat
recon=ifft2c(z3);

a=recon(mask>0);
b=reshape(a,[252,240]);

fb=fft(b,[],2);
ref=load('refg_0_10_30_6cy_TR2.0');
[cc,p]=corrcoef(ref',abs(b(150,:)));

%%

recon=ts_detrend(recon,1:size(recon,4),0);

tmp=reshape(recon,[64*64*9,240]);

a=tmp'*tmp;
[u,d]=eig(a);

u*a