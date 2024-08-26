fid_dir = 'c:\labhome\data\MT_BOLD_CBV\053111\tbgems_01.fid';
root = 'c:\labhome\vnmrsys';
z = read_fid(fullfile(fid_dir,'fid'));
z = squeeze(z);
%z = dcCorr2(z);
if readPar(fid_dir,'ni')==0  %no tabc in vnmrj
    tabfile = readPar(fid_dir,'petable');
    z=tabc(z,fullfile(root,'tablib',tabfile(2:end-1)));
end

[x,y]=meshgrid(-64:63,-64:63);
thr=exp(-(x.^2+y.^2)/2/25.^2);

x=-64:63;
thr=exp(-x.^2/2/25.^2);
thr=repmat(thr,[128,1]);
mask=rand(size(z))<thr;

%z=circshift(z,[64,64]);
%mask=circshift(mask,[64,64]);
%thr=circshift(thr,[64,64]);    %do not use them

fz = fft2(z);
fz2 = fftshift(fz);
figure;imagesc(fliplr(abs(fz2)));
colormap(gray);
%%
z2=fftshift(z);
mask2=fftshift(mask);
data=z2.*mask2;
data=data/max(data(:))*5000;
lambda=0.3;
mu=0.3;
c1=0;
%XFM = Wavelet('Daubechies',4,4);
XFM=discct2;
gamma=mu/1000;

res=split_bregman(mask2,data,mu, lambda, gamma,20, 4,XFM,c1);
res=fftshift(res);
figure(100), imshow(abs(res),[]), drawnow
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fz=ifft2c(z);

data=z.*mask;
N = size(data); 	% image Size
DN = size(data); 	% data Size
TVWeight = 0; 	% Weight for TV penalty
xfmWeight = 0.1;	% Weight for Transform L1 penalty
Itnlim = 10;		% Number of iterations

pdf=thr;
%generate Fourier sampling operator
FT = p2DFTxp(mask, N, 1, 2);

% scale data
im_dc = FT'*(data.*mask./pdf);

figure;imagesc(abs(im_dc));
colormap(gray);

scl=max(abs(im_dc(:)));
data = data/scl;
im_dc = im_dc/scl;

%generate transform operator
XFM = Wavelet_ISU([size(data,1),size(data,2)]);	% Wavelet
wvlt=XFM*fz/scl;
tmp=abs(wvlt(:));
tmp=sort(tmp,'descend');
figure(100); plot(tmp,'r');
set(gca,'yscale','log');
set(gca,'xscale','log');


%XFM=1;
% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOPxp;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

%figure(101), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
for n=1:10
	res = fnlCg_xp(res,param);
	im_res = XFM'*res;
	
    figure(101); imshow(abs(im_res),[]), drawnow
    tmp=abs(res(:));
    tmp=sort(tmp,'descend');
    figure(100); hold on;plot(tmp,'b');
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    
end

