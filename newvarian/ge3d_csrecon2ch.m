function ge3d_csrecon2ch(fid_prefix,twoch,slice)
%ge3d_csrecon2ch(fid_prefix,twoch,slice)
% if specify slice, then only that slice of the first array element on channel 1
% will be reconstructed.
% twoch: true for 2 channel and false for 1 channel.


np=readPar(fid_prefix,'np');
nv=readPar(fid_prefix,'nv');
nv2=readPar(fid_prefix,'nv2');
ni2=readPar(fid_prefix,'ni2');


mask=mask_cs2d(fid_prefix);
m2=shiftdim(mask,-1);
m2=repmat(m2,[np/2,1,1]);

arraydim=readPar(fid_prefix,'arraydim');
ne=arraydim/ni2/2;

z=read_fid([fid_prefix,'.fid']);

if twoch
  z=reshape(z,[size(z,1),size(z,2),2*ne,size(z,3)/2/ne]);
  z=permute(z,[1,2,4,3]);
  z2=zeros(np/2,nv,nv2,2*ne);
  m3=repmat(m2,[1,1,1,2*ne]);
  z2(m3>0)=z;
else
  z=reshape(z,[size(z,1),size(z,2),size(z,3)/ne,ne]);
  z2=zeros(np/2,nv,nv2,ne);
  m3=repmat(m2,[1,1,1,ne]);
  z2(m3>0)=z;  
end

fz=ifft(z2,[],1);
fz=ifft(fz,[],2);
fz=ifft(fz,[],3);

fz=circshift(fz,[np/4,nv/2,nv2/2,0]);


if exist('slice','var')
    
 lpe = readPar(fid_prefix,'lpe');%in cm
 lpe2 = readPar(fid_prefix,'lpe2');%in cm

 pro=readPar(fid_prefix,'pro'); %in cm
 ppe=readPar(fid_prefix,'ppe'); %in cm
 ppe2=readPar(fid_prefix,'ppe2'); %in cm
 
 fz=shift_image_fracvox(fz,[0,-ppe/lpe*nv,-ppe2/lpe2*nv2]);
 tmp=abs(squeeze(fz(slice,:,:,1)));
 figure;imagesc(tmp);
 colormap(gray);
 
 
end

z3=circshift(z2,[np/4,0,0,0]);
fz1=ifft(z3,[],1);
fz1=circshift(fz1,[np/4,0,0,0]);

%fz1=flipdim(fz1,2);   % needed for uFT
%fz1=flipdim(fz1,3);   % needed for uFT

im_cs=zeros(size(fz1));

for islice=1:size(fz1,1)
    
    for nt=1:size(fz1,4)
        if exist('slice','var') 
            if islice == slice && nt==1 
               run_cs( squeeze(fz1(islice,:,:,nt)),mask,false);
               return;
            else
               continue;
            end
        else
          im_cs(islice,:,:,nt)=run_cs( squeeze(fz1(islice,:,:,nt)),mask,true);
        end
   end
end

if ~exist('slice','var')
 
 lro = readPar(fid_prefix,'lro');%in cm
 lpe = readPar(fid_prefix,'lpe');%in cm
 lpe2 = readPar(fid_prefix,'lpe2');%in cm

    pro=readPar(fid_prefix,'pro'); %in cm
    ppe=readPar(fid_prefix,'ppe'); %in cm
    ppe2=readPar(fid_prefix,'ppe2'); %in cm
    sz=size(z);
    delta = [lro*10/sz(1),lpe*10/sz(2),lpe2*10/sz(3)];
        
    orig_delta(1,:) = [pro,ppe,ppe2]*10-delta.*(sz(1:3)/2-[0.5,0.5,0.5]);
    orig_delta(2,:) = delta;
    
    
    orig_delta(:,[1,3])=-orig_delta(:,[1,3]);

    if strcmp(orient,'sag');
       im_cs=permute(im_cs,[3,2,1,4]);
       fz=permute(fz,[3,2,1,4]);
       orig_delta=orig_delta(:,[3,2,1]);
    end

    
    info.ORIGIN = orig_delta(1,:);
    info.DELTA = orig_delta(2,:);
    
    
    fz=shift_image_fracvox(fz,[0,-ppe/lpe*nv,-ppe2/lpe2*nv2]);
    im_cs=shift_image_fracvox(im_cs,[0,-ppe/lpe*nv,-ppe2/lpe2*nv2]);
 
    write_afni(abs(im_cs),info,[fid_prefix,'_cs']);
    write_afni(abs(fz),info,[fid_prefix,'_uFT']);
    
    
end

function im_res=run_cs(z,mask,nodisplay)


data=z.*mask;
data=data/max(abs(data(:)));
N = size(data); 	% image Size
TVWeight = 0.002; 	% Weight for TV penalty
xfmWeight = 0.003;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

%generate Fourier sampling operator
FT = p2DFT(mask, N, 1, 2);

% scale data
im_dc = FT'*data;


data = data/max(abs(im_dc(:)));
im_dc = im_dc/max(abs(im_dc(:)));

%generate transform operator
XFM = Wavelet_rect('Daubechies',4,4);	% Wavelet

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

if ~nodisplay
  figure(100), imagesc(abs(im_dc));colormap(gray);drawnow;
end

res = XFM*im_dc;

% do iterations
for n=1:40
	res = fnlCg_xp(res,param);
	im_res = XFM'*res;
    if n==2
        ph=angle(im_res);
        param.FT=p2DFT(mask,N,cos(ph)+1i*sin(ph));
        im_res=im_res.*(cos(ph)-1i*sin(ph));
        res=XFM*im_res;
    end
    
    if ~nodisplay
	  figure(100), imagesc(abs(im_res));colormap(gray);drawnow;
    end
end



