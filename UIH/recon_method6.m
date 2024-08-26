function [img,img_sos]=recon_method6(fd,motionPar,shotIndex,weight,FOV,ro_pe_par,ds)
%fd should be in the k-space for all dimensions

mask=squeeze(sos(sos(fd,1),4))>1;
[rows,cols]=fullySampledRegion(mask);

mtmp = 0*mask;
mtmp(rows,cols)=1;
mtmp = logical(mtmp(:));

[k,data]=get_k_data_moco_shotIndex_zbb(fd,motionPar,shotIndex,FOV,ro_pe_par,ds);

sel=sos(sos(data,3),1)>0;

weightv=weight(:);

weight2=weightv(sel);

mtmp2=mtmp(sel);
%%
data2=data(:,sel,:);
k2=k(:,sel,:);

%D=reshape(data2,[size(data2,1)*size(data2,2),size(data2,3)]);
%[U,S,V] = svd(D,'econ');
%nCoil = max(find(diag(S)/S(1)>0.05));
%data3 = reshape(D*V(:,1:nCoil),size(data2,1),size(data2,2),nCoil);


ical=size(data,1)/2+1;

%%
calib=(data2(:,mtmp2,:));
kc=(k2(:,mtmp2,:));
weightc=squeeze(weightv(mtmp));
%%
sz=size(fd);
img=zeros(sz);

par.ksize=[7,7]; %bright diagonal region in kro=0 plane when ksize=[6,6];

par.nIterCG=15;

par.imSize = [sz(2),sz(3)];
par.calSize = [10,10];
par.lambda=1;


kc_c = squeeze(kc(ical,:,2:3));
dc_c = reshape(calib(ical,:,:),[size(calib,2),1,size(calib,3)]);
nCoil=size(calib,3);
%
GFFT1 = NUFFT(kc_c,weightc, [0,0] , par.imSize);
im_c = GFFT1'*dc_c;
kData = fft2c(im_c);


kCalib = crop(kData,[par.calSize,nCoil]); % crop center k-space for calibration
CalibTyk = 0.02;
kernel = calibSPIRiT(kCalib, par.ksize, nCoil, CalibTyk);
GOP = SPIRiT(kernel, 'image',par.imSize); % init the SPIRiT Operator
disp('Done Calibrating');
ltime=tic;

sz_data2=size(data2);
%parpool(4);  %parfor produced an error
%{
警告: 类 'SPIRiT' 的元素与当前构造函数定义不匹配。这些元素已转换为结构体。
> 位置：parallel.internal.pool.deserialize (第 33 行)
位置: parallel.internal.pool.deserializeFunction (第 17 行)
位置: remoteParallelFunction (第 29 行)
错误使用 iterapp
用户提供的 function ==> @(x,tflag)afun(x,NUFFTOP,GOP,dataSize,imSize,lambda,tflag) 失败，并返回以下错误:

 输入参数的数目不足。

出错 lsqr (第 180 行)
    u = u - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:},'notransp');

出错 cgNUSPIRiT (第 27 行)
[res,FLAG,RELRES,ITER,RESVEC,LSVEC] = lsqr(@(x,tflag)afun(x,NUFFTOP,GOP,dataSize, imSize,lambda,tflag), b, [], nIter,speye(N,N),speye(N,N), x0(:));
%}
for i=1:size(k2,1)

%for i=1:sz(1)
% k space coordinates
ki = squeeze(k2(i,:,2:3));
% k space data
kdi = reshape(data2(i,:,:),[sz_data2(2),1,nCoil]);

GFFT_u = NUFFT(ki,weight2, [0,0], par.imSize);

%  [~,~,img(i,:,:,:)] = do_SAKE_ESPIRiT_NUFFT(ki,kdi,kci,kcdi,par,ESP);

im_dc = GFFT_u'*kdi;
tmp = cgNUSPIRiT(kdi,im_dc,GFFT_u,GOP,par.nIterCG,par.lambda);
% tmp = cgNUSPIRiT(kdi,im_dc,GFFT_u,GOP,par.nIterCG*4,par.lambda);
%%
img(i,:,:,:)=tmp;
time_left(i,sz(1),toc(ltime));
%         [res,resw,img(i,:,:,:)] = do_SAKE_ESPIRiT_NUFFT(ki,kdata,kci,kcdata,par);
end
img = ifft1c(img,1);

img_sos = sos(img,4);


