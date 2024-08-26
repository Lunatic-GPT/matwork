function [res,resESPIRiT,weights,reskESPIRiT]=do_SAKE_ESPIRiT(DATA,calib,ksize,wnthresh,eigThresh_im,nIterCG )
%DATA: Nx*Ny*Nch;
%calib: calibration data; sx*sy*Nch
%ksize: ESPIRiT kernel-window-size
%wnthresh: Window-normalized number of singular values to threshold
%eigThresh_im: threshold of eigenvectors in image space
%nIterCG: interation number for cgESPIRiT
%% for debug/test
%{
load 8_fse_mx_3d_fatnav2_ACS5_fd170.mat
DATA=fd_170;
[row,col]=fullySampledRegion(sos(DATA,3)>0);
row([156:164,190:198])=1;
col([115:123,149:157])=1;
calib=DATA(row>0,col>0,:);
ksize=[5,5];
wnthresh = 1.8;
 eigThresh_im = 0.9;
 nIterCG=15;
%}
%%
if ~exist('ksize','var')
ksize = [6,6]; 
end

if ~exist('wnthrsh','var')
 wnthresh = 1.8; % Window-normalized number of singular values to threshold
end

if ~exist('eigThresh_im','var')
 eigThresh_im = 0.9; % threshold of eigenvectors in image space
end
if ~exist('nIterCG','var')
nIterCG=15;
end
%%
[sx,sy,Nc] = size(DATA);

sakeIter = 100;
tic; calib = SAKE(calib, ksize, wnthresh,sakeIter, 0);toc
%%
[k,S] = dat2Kernel(calib,ksize);
[M,W] = kernelEig(k(:,:,:,1:floor(wnthresh*prod(ksize))),[sx,sy]);

maps = M(:,:,:,end-1:end);

weights = W(:,:,end-1:end);
weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end-1:end) > eigThresh_im);
weights = -cos(pi*weights)/2 + 1/2;


ESP = ESPIRiT(maps,weights);


DATA=double(DATA);
tic; [reskESPIRiT, resESPIRiT] = cgESPIRiT(DATA,ESP, nIterCG, 0.01,DATA); toc;


res=sos(resESPIRiT.*weights);
