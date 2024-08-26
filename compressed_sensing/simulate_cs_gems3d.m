
fid_prefix='/data/xiaopeng/PhaseImaging_120730/ge3d_08';
z=read_fid([fid_prefix,'.fid']);
fz=fft(z,[],1);
fz=fft(fz,[],2);
fz=fft(fz,[],3);

fz=circshift(fz,[128,128,128]);

figure;imagesc(abs(fz(:,:,128)));
colormap(gray);


%%

n1=256;
n2=64;
xv=linspace(-1/2,1/2-1/n1,n1);
yv=linspace(-1/2,1/2-1/n2,n2);

[x,y]=meshgrid(yv,xv);
thr=exp(-(x.^2+y.^2)/2/0.125.^2);


mask=rand(n1,n2)<thr;
%m2=shiftdim(mask,-1);
%m2=repmat(m2,[256,1,1]);

disp(length(find(mask>0))/n1/n2);


disp(length(find(mask>0)));
petable=fopen('petable_sig0p125','w');
fprintf(petable,'t1 =\n');
ind=find(mask>0);

for i=1:length(find(mask>0))
    fprintf(petable,'  %d\n',ind(i));
end
fclose(petable);
figure;imagesc(mask);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slice=128;

data=squeeze(fz1(128,:,:)).*mask;

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
param.XFM = TVOP;
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
