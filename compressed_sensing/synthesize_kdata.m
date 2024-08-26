function [z3,mask]=synthesize_kdata(nois_scl,bold_scl)
tmp=load('data/kdata_gems_cs_14.mat');
z=mean(tmp.z,4);

recon=ifft2c(z);
%figure;subplot(1,2,1);imshow(abs(recon(:,:,4)),[]);
%subplot(1,2,2);imshow(angle(recon(:,:,4)),[]);

z=z(:,:,4);

%ref=load('refg_0_10_30_6cy_TR2.0');
ref=sin((0:239)/40*2*pi)+1;
ref=ref/max(ref)*bold_scl;
ref=reshape(ref,[1,1,1,length(ref)]);

nt=length(ref);

z=repmat(z,[1,1,1,nt]);

sz=size(z);
z=z./max(abs(z(:)));
ref=repmat(ref,[sz(1:3),1]);

nois=nois_scl*randn_white(length(z(:)));
nois=reshape(nois,size(z));

z2=z+nois;

mask=zeros(64,64,1);
mask(48:52,27:31)=1;
recon2=ifft2c(z2);


mask=repmat(mask,[1,1,1,sz(4)]);
% add BOLD effect;
recon3=recon2.*(1+mask.*ref);

z3=fft2c(recon3);


%{
im=abs(recon2(:,:,1,1));
im=im/max(im(:))*100;
cm_u=gray(100);
cm_o=[1,0,0];
[out,cm_out]=combine_over_under(im,mask,cm_u,cm_o,mask>0);
 figure;imshow(out,cm_out);

figure;subplot(1,2,1);imshow(abs(recon2(:,:,1)),[]);
subplot(1,2,2);imshow(angle(recon2(:,:,1)),[]);

%}


