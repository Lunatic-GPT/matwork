load brain512


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
mask=rand(size(z))<thr;

fz = fft2(z);
fz2 = fftshift(fz);
figure;imagesc(fliplr(abs(fz2)));
colormap(gray);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=z.*mask;
N = size(data); 	% image Size
DN = size(data); 	% data Size
TVWeight = 0.002; 	% Weight for TV penalty
xfmWeight = 0.005;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

pdf=thr;
%generate Fourier sampling operator
FT = p2DFT(mask, N, 1, 2);

% scale data
im_dc = FT'*(data.*mask./pdf);

figure;imagesc(abs(im_dc));
colormap(gray);

data = data/max(abs(im_dc(:)));
im_dc = im_dc/max(abs(im_dc(:)));

%generate transform operator
XFM = Wavelet('Daubechies',4,4);	% Wavelet

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
tic
for n=1:6
	res = fnlCg(res,param);
	im_res = XFM'*res;
	figure(100), imshow(abs(im_res),[]), drawnow
end
toc
