z=PhantomCircle(128,2);
z=repmat(z,[1,1,1,100]);
fz=ifft2c(z);
%fz=fz.*repmat(exp(-1i*pi*(1:128)'/64),[1,128]);
z=fft2c(fz);
z=z/max(abs(z(:)));
nois=randn_white(128*128*100);
nois=reshape(nois,[128,128,1,100]);
z=z+max(abs(z(:)))*0.001*nois;
%m=mask_cs1d('data/gems_cs_03');
im=ifft2c(z);
figure;imshow(abs(im(:,:,1,1)),[]);
a=loadtable('petable_R4_N128_nTR240_ns1');

mask=zeros(128,128,1,100);
for i=1:100
mask(:,a((i-1)*32+1:i*32)+65,1,i)=1;
end

N = size(z); 	% image Size
%generate Fourier sampling operator
FT = p2DFTxp(mask, N, 1, 2);
data=z.*mask;
TVWeight = 0.001; 	% Weight for TV penalty
xfmWeight = 0.001;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

im_dc = FT'*data;
scl=max(abs(im_dc(:)));
data = data/scl;
im_dc = im_dc/scl;

XFM=FTt;%Wavelet_ISU(sz(1:2));

mu=100;
lambda=0;
gamma=mu/1000;

%u = split_bregman(mask,data, mu, lambda, gamma, 10, 15,XFM,mu);
% scale data
FTtmp=p2DFTxp(ones(size(z)), N, 1, 2);
%im_dc=FTtmp'*z;

%im_dc=ones(size(data));
%im_dc=zeros(size(im_dc));
%generate transform operator
%XFM = Wavelet_rect('Daubechies',4,4);	% Wavelet
sz=size(data);
% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOPxp;
param.data = data;
param.TVWeight =0;     % TV penalty 
param.xfmWeight = 0.01;  % L1 wavelet penalty
param.Itnlim = Itnlim;

isl=1;

 % figure(100);subplot(1,2,1); imshow(abs(im_dc(:,:,isl)),[]);
 % subplot(1,2,2); imshow(angle(im_dc(:,:,isl)),[]);

res = XFM*im_dc;
%res=XFM*im;
%save u u
%load u;
%ph=angle(u);  
%{
mu=10;
lambda=mu/10;
gamma=mu/1000;
u = split_bregman(mask,data, mu, lambda, gamma, 2, 5,XFM,0.01,exp(-1i*ph));
%}
% do iterations
for n=1:20
	res = fnlCg_xp(res,param);
    
	im_res = XFM'*res;
    if n==40
        ph=angle(im_res);
        param.FT=p2DFTxp(mask,N,cos(ph)+1i*sin(ph));
        im_res=im_res.*(cos(ph)-1i*sin(ph));
        res=XFM*im_res;
    end
    
	  figure(100);subplot(1,2,1); imshow(abs(im_res(:,:,isl)),[]);
      subplot(1,2,2);imshow(angle(im_res(:,:,isl)),[]);drawnow;
    
end

