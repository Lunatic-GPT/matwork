function ge3d_csrecon(fid_prefix,slice)


tabfile=readPar(fid_prefix,'petable');

if ~exist('tname','var')
    tname = 't1';
end

if strcmp(computer,'PCWIN64')
   fid=fopen(['z:/tmp.home/xiaopeng/newvnmrsys/tablib/',tabfile(2:end-1)]);
else
   fid=fopen(['/home/xiaopeng/vnmrsys/tablib/',tabfile(2:end-1)]);
end
while 1
    a=fgetl(fid);
    if feof(fid)
        error(['Table ',tname,' not found/']);
    end
    if ~isempty(findstr(a,tname))
        tab=[];
        while 1
            b = fgetl(fid);
            if isempty(str2double(b)) || isnan(str2double(b))
                break;
            end
            tab(end+1)=str2double(b);
            if  feof(fid)
                break;
            end
        end
        break;
    end
    
end
cs_pe=readPar(fid_prefix,'cs_pe');
cs_pe2=readPar(fid_prefix,'cs_pe2');
np=readPar(fid_prefix,'np');
nv=readPar(fid_prefix,'nv');
nv2=readPar(fid_prefix,'nv2');

tab=tab(1:cs_pe*cs_pe2);
mask=zeros(nv,nv2);
mask(tab)=1;
m2=shiftdim(mask,-1);
m2=repmat(m2,[np/2,1,1]);

z=read_fid([fid_prefix,'.fid']);
z=z(:,:,1:2:end-1);
z2=zeros(np/2,nv,nv2);
z2(m2>0)=z;

fz=fft(z2,[],1);
fz=fft(fz,[],2);
fz=fft(fz,[],3);

fz=circshift(fz,[np/4,nv/2,nv2/2]);

figure;imagesc(abs(squeeze(fz(slice,:,:))));
colormap(gray);

z3=circshift(z2,[np/4,0,0]);
fz1=fft(z3,[],1);
fz1=circshift(fz1,[np/4,0,0]);



im_cs=zeros(size(fz1));

figure;imagesc(abs(im_dc));
colormap(gray);


for slice=1:size(fz1,1)
    fprintf('Slice %d\n',slice);
    
data=squeeze(fz1(slice,:,:)).*mask;

N = size(data); 	% image Size
TVWeight = 0.002; 	% Weight for TV penalty
xfmWeight = 0.005;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

xv=linspace(-1/2,1/2-1/nv,nv);
yv=linspace(-1/2,1/2-1/nv2,nv2);

[x,y]=meshgrid(yv,xv);
%generate Fourier sampling operator
FT = p2DFT(mask, N, 1, 2);

% scale data
im_dc = FT'*(data.*mask);


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

figure(100), imagesc(abs(im_dc));colormap(gray);

res = XFM*im_dc;

% do iterations
tic
for n=1:4
	res = fnlCg(res,param);
	im_res = param.XFM'*res;
	figure(100), imagesc(abs(im_res));colormap(gray);
end
toc;
end

