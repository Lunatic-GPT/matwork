function [z3,mask]=synthesize_kdata_mask2_nih(nois_scl,bold_scl)
% synthesize_kdata_mask2_nih(nois_scl,bold_scl)

tmp=load('kdata_gems_cs_14.mat');
z=mean(tmp.z,4);

recon=ifft2c(z);
%figure;subplot(1,2,1);imshow(abs(recon(:,:,4)),[]);
%subplot(1,2,2);imshow(angle(recon(:,:,4)),[]);

%z=z(:,:,4);

ref=load('refg_0_10_30_6cy_TR2.0');
ref=ref/max(ref)*0.02;
%ref=sin((0:239)/40*2*pi)+1;
ref=ref/max(ref)*bold_scl;
nt=length(ref);


mask=zeros(64,64,9,nt);
mask(19:25,45:50,4:6,:)=1;
mask(45:51,43:48,4:6,:)=1;

% add BOLD effect;
sz=size(z);
ref=reshape(ref,[1,1,1,nt]);
ref=repmat(ref,[sz(1:3),1]);
recon=repmat(recon,[1,1,1,nt]);
recon2=recon.*(1+mask.*ref);
nois=nois_scl*randn_white(length(recon2(:))).*exp(1i*2*pi*rand(1,length(recon2(:))));
nois=reshape(nois,size(recon2));

recon3=recon2/max(recon2(:))+nois;

z_true=fft2c(recon3-nois);
z3=fft2c(recon3);


save(sprintf('synthesize_kdata_mask2_nih_nois%4.3f_bold%4.3f.mat',nois_scl,bold_scl),'z3','z_true','mask');

write_afni(abs(recon3),sprintf('synthesize_kdata_nih_nois%4.3f_bold%4.3f',nois_scl,bold_scl));

%write_afni(mask,'act_mask2');
%%
%{
ud=abs(recon2(:,:,4));
ud=ud/max(ud(:))*100;
ov=mask(:,:,4);

[out,cm_out]=combine_over_under(ud,ov,gray(100),[1,0,0],ov>0);
figure;imshow(out,cm_out);
%}


