function recon_gmems(a,format)
%gems_csrecon(a,reduction,format)

if ~exist('format','var')
    format='a';
end
tstart=tic;

of1=readbPar(fullfile(a,'method'),'PVM_SPackArrPhase1Offset');
of0=readbPar(fullfile(a,'method'),'PVM_SPackArrReadOffset');
fov=readbPar(fullfile(a,'method'),'PVM_Fov');

ns=readbPar(fullfile(a,'method'),'PVM_SPackArrNSlices');

data=readb_fid(a);

pe=readbPar(fullfile(a,'acqp'),'ACQ_spatial_phase_1');
NR = readbPar(fullfile(a,'method'),'PVM_NRepetitions');


mtrx=readbPar(fullfile(a,'method'),'PVM_Matrix');
pe=round(pe*mtrx(2)/2)+mtrx(2)/2+1;
ne=readbPar(fullfile(a,'method'),'PVM_NEchoImages');

obj=readbPar(fullfile(a,'method'),'PVM_ObjOrderList');

[nch,ch_scl]=nChan(a);


iNR = length(pe)/mtrx(2);
%data=reshape(data,[mtrx(1),nch,ns,mtrx(2)/reduction,reduction*iNR,NR/iNR]);
data=reshape(data,[mtrx(1),nch,ne,ns,mtrx(2),NR]);

data(:,:,2:2:end,:,:,:)=flipdim(data(:,:,2:2:end,:,:,:),1);
data(:,:,:,obj+1,:,:)=data;



data=permute(data,[1,5,4,3,6,2]);


img=ifft2c(data);

    fprintf('gems_csrecon finished in %4.1f s\n',toc(tstart));

    
    
function d2=phaseCorr(d2,nav)

  ph=angle(nav);
  
  ph=repmat(ph,[size(d2,1),1,1,1,1]);
  
 d2=d2.*exp(-1i*ph);

 
    function out=omp_mat(b,m)
        
        out=fft2c(ifft1c(b,4)).*m; 
        out=out(:);
        
        
   function out=omp_matt(b,m)
            
        d=zeros(size(m));
        d(m>0)=b;
        out=ifft2c(fft1c(d,4));
        out=out(:);        
        
        
function im_res=run_cs_FTt(z,mask,nodisplay,isl,xfm,niter,xfmWeight,TVWeight)


if ~exist('niter','var')
niter=2;
end
if ~exist('xfmWeight','var')
xfmWeight=0.1;
end
if ~exist('TVWeight','var')
 TVWeight=0.01;
end
z=z/max(abs(z(:)));
data=z.*mask;
%data=data/max(abs(data(:)));
N = size(data); 	% image Size
Itnlim = 8;		% Number of iterations

%generate Fourier sampling operator
FT = p2DFTxp(mask, N, 1, 2);

% scale data

%data = data/max(abs(im_dc(:)));
%im_dc = im_dc/max(abs(im_dc(:)));

%generate transform operator
%XFM = Wavelet_rect('Daubechies',4,4);	% Wavelet
%sz=size(data);
%XFM=Wavelet_ISU(sz(1:2));
% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.data = data;
param.TVWeight =TVWeight;%TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;%xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;


in=sum(data,4);
nm=sum(mask,4);
nm(nm==0)=1;
in=in./nm;
im_dc=FT'*(repmat(in,[1,1,1,size(data,4)]));
%im_dc=ifft2c(repmat(in,[1,1,1,size(data,4)]));

sd=std(data,0,4);
sd2=sd(nm==size(mask,4));

fprintf('average noise = %f\n',sqrt(mean(sd2.^2)));
%figure;hist(sd2,0:0.002:0.1);drawnow;


if strmatch(xfm,'ft')
XFM=FTt;
elseif strmatch(xfm,'klt')
 m2=(nm==size(mask,4));
 d3=zeros(length(find(m2>0)),size(mask,4));
for i=1:size(mask,4)
    tmp=data(:,:,:,i);
    d3(:,i)=tmp(m2);
end
 XFM=KLT(d3);
end

sz=size(data);
param.TV=Wavelet_ISU(sz(1:2));

param.XFM=XFM;

if ~nodisplay
  figure(100);
  for i=1:length(isl)
    subplot(1,length(isl),i); imshow(abs(im_dc(:,:,isl(i))),[]); 
    drawnow;   
  end
      
 end

res = param.XFM*im_dc;
% do iterations
for n=1:niter
	res = fnlCg_xp(res,param);
	im_res = param.XFM'*res;
    if n==40
        ph=angle(im_res);
        param.FT=p2DFTxp(mask,N,cos(ph)+1i*sin(ph));
        im_res=im_res.*(cos(ph)-1i*sin(ph));
        res=param.XFM*im_res;
    end
    
    if ~nodisplay
	  figure(100);
      for i=1:length(isl)
         subplot(1,length(isl),i); imshow(abs(im_res(:,:,isl(i))),[]);
         drawnow;   
      end
      %subplot(1,2,2);imshow(angle(im_res(:,:,isl)),[]); drawnow;
    end
end

function im_res=run_cs_1TR(z,mask,nodisplay,isl)


data=z.*mask;
data=data/max(abs(data(:)));
N = size(data); 	% image Size
TVWeight = 0.003; 	% Weight for TV penalty
xfmWeight = 1;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

%generate Fourier sampling operator
FT = p2DFTxp(mask, N, 1, 2);

% scale data
im_dc = FT'*data;


data = data/max(abs(im_dc(:)));
im_dc = im_dc/max(abs(im_dc(:)));

%generate transform operator
%XFM = Wavelet_rect('Daubechies',4,4);	% Wavelet
sz=size(data);
XFM=Wavelet_ISU(sz(1:2));
% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOPxp;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;


if ~nodisplay
  figure(100);subplot(1,2,1); imshow(abs(im_dc(:,:,isl)),[]);
  subplot(1,2,2); imshow(angle(im_dc(:,:,isl)),[]);
end

res = XFM*im_dc;

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
    
    if ~nodisplay
	  figure(100);subplot(1,2,1); imshow(abs(im_res(:,:,isl)),[]);
      subplot(1,2,2);imshow(angle(im_res(:,:,isl)),[]);
    end
end

       