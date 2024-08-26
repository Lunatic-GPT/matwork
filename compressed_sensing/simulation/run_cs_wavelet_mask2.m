function run_cs_wavelet_mask2(nois,tvWeight,ftWeight,m,suffix)

bold=0.015;
tmpname = sprintf('synthesize_kdata_mask2_nih_nois%4.3f_bold%4.3f.mat',nois,bold);
if ~exist(tmpname,'file')
[z2,mask]=synthesize_kdata_mask2(nois,0.02);
else
tmp=load(tmpname);
z2=tmp.z3;
mask=tmp.mask;
end

%z2=tmp.z2;
%recon_orig=ifft2c(z2);

z2=m.*z2;

recon=run_cs(z2,m,true,4:6,tvWeight,ftWeight);        
        
        prefix = sprintf('run_cs_wavelet_nih_nois%4.3f_bold%4.3f_tv%4.3f_ft%4.3f%s',nois,bold,tvWeight,ftWeight,suffix);
        write_afni(abs(recon),prefix);
        write_afni(angle(recon),[prefix,'_ph']);
        
  % recon=ifft2c(z2);
        %mask=zeros(64,64);
        %mask(48:52,27:31)=1;
      %  mask(50,29)=1;
%        ts=mean_roi(abs(recon),mask(:,:,:,1));
%        ts_orig=mean_roi(abs(recon_orig),mask(:,:,:,1));
%{        
        figure;subplot(2,1,1);plot(ts);
        
        subplot(2,1,2);plot(ts_orig);
        
        figure;subplot(1,2,1);imshow(abs(recon(:,:,6,1)),[]);
        subplot(1,2,2);imshow(abs(recon_orig(:,:,6,1)),[]);
 %}       
        
function im_res=run_cs(z,mask,nodisplay,isl,tvWeight,ftWeight)

data=z.*mask;
data=data/max(abs(data(:)));
N = size(data); 	% image Size


%generate Fourier sampling operator


%generate transform operator
%XFM = Wavelet_rect('Daubechies',4,4);	% Wavelet
sz=size(data);
XFM=Wavelet_ISU(sz(1:2));
%XFM=FTt;
% initialize Parameters for reconstruction
param = init;
param.XFM = XFM;
param.TV = TVOPxp;
param.TVWeight =tvWeight;     % TV penalty 
param.xfmWeight = ftWeight;  % L1 wavelet penalty
param.Itnlim = 8;


in=sum(data,4);
nm=sum(mask,4);
nm(nm==0)=1;
in=in./nm;
im_dc=ifft2c(repmat(in,[1,1,1,size(data,4)]));

%tmp=load('ktFOCUSS_default_sin0p04');
%im_dc=tmp.recon;
%tmp=load('recon_orig_sin0p04');
%im_dc=tmp.recon_orig;
%param.data=param.data;
%im_dc=im_dc*max(abs(im_dc_orig(:)))/max(abs(im_dc(:)));
%scale=max(abs(im_dc(:)));
%im_dc=im_dc/scale;
%param.data=param.data/scale;
ns=length(isl);
if ~nodisplay
  figure(100);
  for i=1:ns
   subplot(1,ns,i); imshow(abs(im_dc(:,:,i)),[]);drawnow;
  end
end

res = XFM*im_dc;
% do iterations

for it=1:size(res,4)
    fprintf('Recon image %d\n',it);
    
param.data = data(:,:,:,it);
FT = p2DFTxp(mask(:,:,:,it), N, 1, 2);
param.FT=FT;
for n=1:8
	res(:,:,:,it) = fnlCg_xp(res(:,:,:,it),param);
	im_res = XFM'*res(:,:,:,it);
    if n==40
        ph=angle(im_res);
        param.FT=p2DFTxp(mask,N,cos(ph)+1i*sin(ph));
        im_res=im_res.*(cos(ph)-1i*sin(ph));
        res=XFM*im_res;
    end
    
    if ~nodisplay
	   figure(100);
     for i=1:ns
      subplot(1,ns,i); imshow(abs(im_res(:,:,i)),[]);drawnow;
     end
    end
end

end

im_res=XFM'*res;