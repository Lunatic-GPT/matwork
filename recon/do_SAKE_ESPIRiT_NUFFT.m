function [res,resw,resMultiCoil]=do_SAKE_ESPIRiT_NUFFT(k,Data,kc,Datac,par)
%Data: N*Nch;
%k: N*2

%par contains the following fields:
%ksize,wnthresh,eigThresh_im,nIterCG,cind
%%

ksize = par.ksize; 
wnthresh = par.wnthresh; % Window-normalized number of singular values to threshold
eigThresh_im = par.eigThresh_im; % threshold of eigenvectors in image space
nIterCG=par.nIterCG;
nIterSAKE=par.nIterSAKE;
imSize=par.imSize;
calSize=par.calSize;

%{
load 8_fse_mx_3d_fatnav2_ACS5_fd170.mat;
DATA=reshape(fd_170,[352,270,1,64]);
[row,col]=fullySampledRegion(sos(DATA,4)>0);
row([156:164,190:198])=1;
col([115:123,149:157])=1;


ksize=[6,6];
wnthresh = 1.8;
eigThresh_im = 0.9;
for nIterCG=[15,30]
nIterSAKE=100;
imSize=[352,270];
calSize=[sum(row),sum(col)];

[k,Data]=get_k_data_moco(DATA,[],[200,200,2],[1,2,3],[0,0,0]);
sel=sos(Data,2)>0;
k=k(sel,:);
Data=Data(sel,:);

[kc,Datac]=get_k_data_moco(DATA(row,col,:,:),[],[200,200,2],[1,2,3],[0,0,0]);
sel=sos(Datac,2)>0;
kc=kc(sel,:);
Datac=Datac(sel,:);
%}
%%
%{
load brain_8ch.mat;

ncalib = 48;
mask = mask_nocalib_x2half;     % choose a 3x no calibration mask
DATAc = DATA.* repmat(mask,[1,1,size(DATA,3)]);
calibc = crop(DATAc,[ncalib,ncalib,8]);

DATAc=reshape(DATAc,[200,200,1,8]);
calibc=reshape(calibc,[ncalib,ncalib,1,8]);

ksize=[6,6];
wnthresh = 1.8;
eigThresh_im = 0.9;
nIterCG=15;
nIterSAKE=100;

imSize=[200,200];
calSize=[ncalib,ncalib];

[k,Data]=get_k_data_moco(DATAc,[],[200,200,2],[1,2,3],[0,0,0]);
sel=sos(Data,2)>0;
k=k(sel,:);
Data=Data(sel,:);

[kc,Datac]=get_k_data_moco(calibc,[],[200,200,2],[1,2,3],[0,0,0]);
sel=sos(Datac,2)>0;
kc=kc(sel,:);
Datac=Datac(sel,:);
%}
%%

calib = SAKE_nufft(Datac, kc(:,1:2), calSize,ksize, wnthresh,nIterSAKE);

%%

knl = dat2Kernel(calib,ksize);

[M,W] = kernelEig(knl(:,:,:,1:floor(wnthresh*prod(ksize))),imSize);


nmaps=2;
maps = M(:,:,:,end-nmaps+1:end);

weights = W(:,:,end-nmaps+1:end);
weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end-nmaps+1:end) > eigThresh_im);
weights = -cos(pi*weights)/2 + 1/2;

ESP = ESPIRiT(maps,weights);



nufft=NUFFT(k(:,1:2),1,[0,0],imSize);

x0=ESP'*(nufft'*Data);

Data2=reshape(Data,[size(Data,1),1,size(Data,2)]);

res = cgL1ESPIRiT(Data2, x0, nufft, ESP, nIterCG, 1, 0, 0,1); % return reconstructed images


resMultiCoil = ESP*res;
%%

 resw=sos(res.*weights,3);

%  save(sprintf('resw_%diter.mat',nIterCG),'resw');
end







