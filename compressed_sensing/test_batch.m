function test_batch

tmp=load('synthesize_kdata_nois0p005_BOLD0p04_sin');
z2=tmp.z2;
recon_orig=ifft2c(z2);
a=loadtable('petable_R4_N64_nTR240_ns1');
m=zeros(64,64,1,240);
for i=1:240
    m(:,a((i-1)*16+1:i*16)+33,1,i)=1;
end

z2=m.*z2;
recon=zeros(size(z2));
        for i=1:size(z2,3)            
          fprintf('Slice %d\n',i); 
          recon(:,:,i,:)=run_cs(z2,m,false,1);        
        end
        
        write_afni(abs(recon),'test_batch');
  % recon=ifft2c(z2);
        mask=zeros(64,64);
        mask(48:52,27:31)=1;
      %  mask(50,29)=1;
        ts=mean_roi(abs(recon),mask);
        ts_orig=mean_roi(abs(recon_orig),mask);
        
        figure;subplot(2,1,1);plot(ts);
        
        subplot(2,1,2);plot(ts_orig);
        
        figure;subplot(1,2,1);imshow(abs(recon(:,:,1,1)),[]);
        subplot(1,2,2);imshow(abs(recon_orig(:,:,1,1)),[]);
        
        
function im_res=run_cs(z,mask,nodisplay,isl)

data=z.*mask;
%data=data/max(abs(data(:)));
N = size(data); 	% image Size
TVWeight = 0; 	% Weight for TV penalty
xfmWeight = 0.001;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

%generate Fourier sampling operator
FT = p2DFTxp(mask, N, 1, 2);

% scale data
im_dc = FT'*data;


%data = data/max(abs(im_dc(:)));
%im_dc = im_dc/max(abs(im_dc(:)));

%generate transform operator
%XFM = Wavelet_rect('Daubechies',4,4);	% Wavelet
%sz=size(data);
%XFM=Wavelet_ISU(sz(1:2));
XFM=FTt;
% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOPxp;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;


in=sum(data,4);
nm=sum(mask,4);
nm(nm==0)=1;
in=in./nm;
im_dc=FT'*repmat(in,[1,1,1,size(data,4)]);

%tmp=load('ktFOCUSS_default_sin0p04');
%im_dc=tmp.recon;
tmp=load('recon_orig_sin0p04');
im_dc=tmp.recon_orig;
%param.data=param.data;
%im_dc=im_dc*max(abs(im_dc_orig(:)))/max(abs(im_dc(:)));
%scale=max(abs(im_dc(:)));
%im_dc=im_dc/scale;
%param.data=param.data/scale;
if ~nodisplay
  figure(100);subplot(1,2,1); imshow(abs(im_dc(:,:,isl)),[]);
  subplot(1,2,2); imshow(angle(im_dc(:,:,isl)),[]); drawnow;
end

res = XFM*im_dc;
% do iterations
for n=1:4
	res = fnlCg_xp(res,param);
	im_res = XFM'*res;
    if n==40
        ph=angle(im_res);
        param.FT=p2DFTxp(mask,N,cos(ph)+1i*sin(ph));
        im_res=im_res.*(cos(ph)-1i*sin(ph));
        res=XFM*im_res;
    end
    
    if ~nodisplay
	  figure(100);subplot(1,2,1); imshow(abs(im_res(:,:,isl)),[]);
      subplot(1,2,2);imshow(angle(im_res(:,:,isl)),[]); drawnow;
    end
end
%im_res=im_res;
        