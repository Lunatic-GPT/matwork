function res=do_ESPIRiT(DATA,calib,ksize,wnthresh,eigThresh_im,nIterCG )
%DATA: Nx*Ny*Nch;
%ksize: ESPIRiT kernel-window-size
%wnthresh: Window-normalized number of singular values to threshold

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

[sx,sy,Nc] = size(DATA);

[k,S] = dat2Kernel(calib,ksize);
[M,W] = kernelEig(k(:,:,:,1:floor(wnthresh*prod(ksize))),[sx,sy]);

maps = M(:,:,:,end-1:end);

weights = W(:,:,end-1:end) ;
weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end-1:end) > eigThresh_im);
weights = -cos(pi*weights)/2 + 1/2;

ESP = ESPIRiT(maps,weights);


DATA=double(DATA);
tic; [reskESPIRiT, resESPIRiT] = cgESPIRiT(DATA,ESP, nIterCG, 0.01,DATA); toc;

img_zf=ifft1c(ifft1c(DATA,1),2);

res=sos(resESPIRiT.*weights);
