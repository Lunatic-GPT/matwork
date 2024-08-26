function epiShaperCS(fid_dir,method,nophase)
% epiShaper32(fid_dir[,method,nophase])
% nophase: do not save phase image.default true.


cs=readPar(fid_dir,'cs');
if ~strcmp(cs(2),'y')
  error('not an CS scan');
end

if ~exist('nophase','var')
    nophase = false;
end

if ~exist(fid_dir,'dir')
  fid_dir=[fid_dir,'.fid'];
end

if ~exist('method','var')
    method='ktf';
end

[pathstr,prefix]=fileparts(fid_dir);

nnav=readPar(fid_dir,'nnav');
nread=readPar(fid_dir,'nread');
nphase=readPar(fid_dir,'nphase');
nseg=readPar(fid_dir,'nseg');
kzero=readPar(fid_dir,'kzero');
ns = readPar(fid_dir,'ns');

if ns>1
    warning('assuming the same mask for all slices');
end
 
im = readPar(fid_dir,'image');
os = readPar(fid_dir,'oversample');

fblipfactor=readPar(fid_dir,'fblipfactor');
fblipfactor=fblipfactor(2:end-1);
if exist(fblipfactor,'file')
    blipfactor=load(fblipfactor);
else
    blipfactor=load(['~/vnmrsys/tablib/',fblipfactor]);
end


%blipfactor=readPar(fid_dir,'blipfactor');


reduction=readPar(fid_dir,'reduction');
nTRcs=readPar(fid_dir,'ncspat');
if isempty(nTRcs)
    nTRcs=1;
end
na=readPar(fid_dir,'arraydim');

mask=zeros(nread/2*os,nphase*reduction,ns,na);
blipfactor=reshape(blipfactor(1:(nphase-1)*nTRcs),[nphase-1,nTRcs]);
blipfactor=cat(1,ones(1,nTRcs),blipfactor);
ind=cumsum(blipfactor,1);
pss=readPar(fid_dir,'pss');
pro=readPar(fid_dir,'pro');
ppe=readPar(fid_dir,'ppe');
lro=readPar(fid_dir,'lro');
lpe=readPar(fid_dir,'lpe');

for j=1:na
    iTR = mod(j-1,nTRcs)+1;
    mask(:,ind(:,iTR),:,j)=1;
end

%save([strtok(fid_dir,'.'),'_mask'],'mask');
%return;
na=length(im);
z_tmp = read_fid(fid_dir);
z=reshape(z_tmp,[nread*os/2,nphase/nseg+nnav+1,ns,nseg,na]);
z=z(:,nnav+2:end,:,:,:);
z=permute(z,[1,2,4,3,5]);
z=reshape(z,[size(z,1),nphase,ns,na]);
z = convertTraces(z);  

z=fov_shift_kspace(z,[pro,0],[lro,lpe]);

z2 = z(:,:,:,im>0);   %implement later
zref=z(:,:,:,im==0);
mask3=mask(:,:,:,im>0);

seqcon=readPar(fid_dir,'seqcon');

if seqcon(3)=='s'
tmp=unique(pss);
ns=length(tmp);
na=size(z2,4)/ns;
z2=reshape(z2,[nread*os/2,nphase,ns,na]);
mask3=reshape(mask3,[nread*os/2,nphase*reduction,ns,na]);

  if size(zref,4)==1  %size(zref,3) is always 1
      zref=repmat(zref,[1,1,ns]);
  else
      zref=reshape(zref,[nread*os/2,nphase,ns]);  % reference scan acquired for each slice.  Assuming the same slice order as the image scans.
  end
end
    
z2= fidPhaseCorr(z2,zref);   %implement later

z3=zeros(nread*os/2,nphase*reduction,ns,size(z2,4));
z3(mask3>0)=z2;

z3=fov_shift_kspace(z3,[0,-ppe],[lro,lpe]);
    
%z2=z(:,:,:,im>0);


%%

switch method
    case 'single'
        
          for j=1:size(z3,4)
              
              %for i=1:size(z3,3)
                img(:,:,:,j)=run_cs(z3(:,:,:,j),mask3(:,:,:,j),false,[1,5,7]);
              %end
          end
    case 'ktf'
        %{
        for j=1:size(z3,3)
            fprintf('Slice %d \n',j);
         tmp=z3(:,:,j,:);
          tmp=permute(tmp,[2,1,3,4]);
          tmp=circshift(tmp,[size(tmp,1)/2,size(tmp,2)/2,0,0]);
          tmp=squeeze(tmp);
          
          recon(:,:,j,:)=cart_ktFOCUSS_KTFOCUSS(tmp);
        end
        
        recon=circshift(recon,[size(z3,1)/2,size(z3,2)/2,0,0]);
         recon=flipdim(recon,1);
        recon=flipdim(recon,2);
        img=permute(recon,[2,1,3,4]);
        %}
        img=zeros(size(mask3));
        for i=1:size(z3,3)            
          fprintf('Slice %d\n',i);
          img(:,:,i,:)=ktFOCUSS_new(z3(:,:,i,:),[],[],3,mask3(:,:,i,:),0.5,0.01,60,2,1);  
        end
        
   case 'ktfklt'
       
        img=zeros(size(z3));
        
        if ~exist([prefix,'_ktf_mag+orig.HEAD'],'file') || ~exist([prefix,'_ktf_ph+orig.HEAD'],'file') 
          prev=zeros(size(z2));
          for i=1:size(z2,3)            
            fprintf('Slice %d\n',i);
            prev(:,:,i,:)=ktFOCUSS_new(z3(:,:,i,:),[],[],3,mask3(:,:,i,:),0.5,0.01,60,2,1);  
   
          end
        else
            mag=BrikLoad([prefix,'_ktf_mag+orig.HEAD']);
            ph=BrikLoad([prefix,'_ktf_ph+orig.HEAD']);
            prev=mag.*exp(1i*ph);
        end
            
        
        for i=1:size(z2,3)            
          fprintf('Slice %d\n',i);
          for iklt=1:3
           img(:,:,i,:)=ktFOCUSS_KLT_new(z3(:,:,i,:),prev,[],[],mask3(:,:,i,:),0.5,0.01,60,2,1);  
           prev=img(:,:,i,:);
          end
        end
           
    case 'ft'
        
        XFMw=0.003;
        TVw=0;
        niter=5;
        for i=1:size(z3,3)
            fprintf('slice %d\n',i);
            img(:,:,i,:)=run_cs_FTt(z3(:,:,i,:),mask3(:,:,i,:),true,1,'ft',niter,XFMw,TVw);
        end
    case 'klt'
         XFMw=0.01;
        TVw=0;
        niter=5;
        for i=1:size(z3,3)
          img(:,:,i,:)=run_cs_FTt(z3(:,:,i,:),mask3(:,:,i,:),true,1,'klt',niter,XFMw,TVw);
        end
        
    case 'modcs'
             
        XFM=Wavelet_ISU([size(z3,1),size(z3,2)]);
        for i=1:size(z3,3)
         img(:,:,i,:)=modcsresCausalDetection_xp(z3(:,:,i,:),mask3(:,:,i,:),XFM,0.05);
        end
    otherwise
        error('Method unknow');
end

% segmentation re-order


pss=readPar(fid_dir,'pss'); %in cm
pss=unique(pss);
b = slice_reorder(img,pss); 

orient=readPar(fid_dir,'orient');
orient=orient(2:end-1);

switch orient 
    case 'trans90'  %1,2,3 in b - l2r,u2d in stimulate - l2r, a2p,s2i in vnmrj
   
   b=flipdim(b,2);  
  % b=flipdim(b,1);  
  order = {'L2R','A2P','S2I'};
   
  case 'trans'
   b = permute(b,[2,1,3,4]);   %do I really need this? yes 
  % b=flipdim(b,2); 
  % b=flipdim(b,1);
   order={'L2R','A2P','S2I'};
   
 case 'sag90'  % 1,2 in b - l2r, u2d in stimulate - a2p, s2i in vnmrj
  %  b = permute(b,[2,1,3,4]);
  
      b=flipdim(b,1);    
  %    b=flipdim(b,2);
      order = {'A2P','S2I','L2R'};
      
    case 'sag' % 1,2 in b - l2r, u2d in stimulate - a2p, s2i in vnmrj
      b=permute(b,[2,1,3,4]);   
      b=flipdim(b,2); 
      b=flipdim(b,1);
      order = {'A2P','S2I','L2R'};
      
    case 'cor90'  
      b=  flipdim(b,1);
   %   b=flipdim(b,2);
      order = {'L2R','S2I','P2A'};
    case 'cor'
      b=permute(b,[2,1,3,4]);
     b=flipdim(b,2);
       b=flipdim(b,1);
      order = {'L2R','S2I','P2A'};
      
    otherwise
        
end

[sgn,pm]=reorient_sdt2afni(order);

tmp=b;
for i=1:3
    if sgn(i)==-1
        tmp=flipdim(tmp,i);
    end
end

afni=permute(tmp,[pm,4]);
    

writesdt4(abs(b),[prefix,'_',method,'_mag']);
if ~nophase
 writesdt4(angle(b),[prefix,'_',method,'_ph']);
end
%%



pro=readPar(fid_dir,'pro');
ppe=readPar(fid_dir,'ppe');
lro=readPar(fid_dir,'lro');
lpe=readPar(fid_dir,'lpe');
nread=nread*os/2;
nphase=nphase*reduction;
lro=lro*os;

thk = readPar(fid_dir,'thk');%in mm
tr = readPar(fid_dir,'tr');
pss=readPar(fid_dir,'pss'); %in cm
%img = slice_reorder(img,pss);    

    sz=size(afni);
    if length(pss)>1
        pss_sort = sort(pss);
        thk = (pss_sort(2)-pss_sort(1))*10;  %pss in cm 
    end
    
    
       
    %varian convention:
    %              axial;  axial90;  coronal; coronal90; sag;  sag90   
    %- to +    pro: A2P;    L2R;      S2I;     L2R;      S2I;  A2P;
    %- to +    ppe: R2L;    A2P;      R2L;     S2I;      P2A;  S2I;
    %- to +    pss: S2I;    S2I;      P2A;     P2A;      L2R;  L2R;
    %  
    %     
    % varian lab coordinates: [x,y,z] = [R2L,P2A,S2I] in roi panel;
    [sgn,pm]=reorient_varian2afni(orient);
    center = [pro,ppe,mean(pss)]*10;
    delta = [lro*10/nread,lpe*10/nphase,thk];
  
    center= center.*sgn;
    center=center(pm);
    
    delta=delta(pm);
    delta([1,3])=-delta([1,3]);
    %    
  %   [img,orig_delta] = reorient_data(img,orig_delta,orient(2:end-1));
  %   should not need to reorient if data is from sdt file.  6/4/2012
    info.ORIGIN = center-delta.*(sz(1:3)/2-[0.5,0.5,0.5]);
    info.DELTA = delta;
    info.TAXIS_FLOATS = [0,tr,0,0,0];
    info.ORIENT_SPECIFIC = [1 3 5]; %L2R,A2P,S2I  % afni convention + is R2L, A2P, I2S
    write_afni(abs(afni),[prefix,'_',method,'_mag'],info);
    if ~nophase
     write_afni(angle(afni),[prefix,'_',method,'_ph'],info);
    end
    
    

%%


function z = convertTraces(z_tmp)
 % reverse the direction of even number traces
 z = zeros(size(z_tmp));
 
 z(:,1:2:end-1,:,:) = z_tmp(:,1:2:end-1,:,:);
 z(1:end,2:2:end,:,:) = z_tmp(end:-1:1,2:2:end,:,:);

function z = dcCorr(z_tmp)
        
 NUM_DC_SUM = 8;
 ave1= mean(z_tmp(1:NUM_DC_SUM/2,:,:,:),1);
 ave2= mean(z_tmp(end-NUM_DC_SUM/2+1:end,:,:,:),1);

 ave=(ave1+ave2)/2;
        
 ave3 = mean(ave,2);       
 z=z_tmp-repmat(ave3,[size(z_tmp,1),size(z_tmp,2),1,1]);
        
 
function z=fidPhaseCorr(z,zref)
  %
        sz = size(z);
        
        fzref = fft(zref,[],1);
        fzref = fzref./abs(fzref);  % normalize
        sz_rep = sz;
        if length(sz)==4  
          sz_rep(1:end-1)=1;
          fzref = repmat(fzref,sz_rep);
        end
        fz = fft(z,[],1);
        fz = fz.*conj(fzref);
        
        z=ifft(fz,[],1);
        z=fftshift(z,1);
        
       

function im_res=run_cs(z,mask,nodisplay,isl)

if ~exist('isl','var')
    isl=1;
end

data=z.*mask;
data=data/max(abs(data(:)));
N = size(data); 	% image Size
TVWeight = 0.05; 	% Weight for TV penalty
xfmWeight = 0;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations

%generate Fourier sampling operator
FT = p2DFTxp(mask, N, 1, 2);

% scale data
im_dc = FT'*data;


data = data/max(abs(im_dc(:)));
im_dc = im_dc/max(abs(im_dc(:)));

%generate transform operator

%if size(z,1)~=size(z,2)
%  XFM = Wavelet_rect('Daubechies',4,4);	% Wavelet
%else

  XFM=Wavelet_ISU([size(z,1),size(z,2)]);
%end

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOPxp;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;
im_dc=zeros(size(im_dc));
if ~nodisplay
  figure(101);
  for i=1:length(isl)
  subplot(length(isl),2,(i-1)*2+1);imagesc(abs(im_dc(:,:,isl(i))));colormap(gray);
  
  subplot(length(isl),2,(i-1)*2+2);imagesc(angle(im_dc(:,:,isl(i))));
  drawnow;
  end
end

res = XFM*im_dc;
%res=XFM*zeros(size(im_dc));
% do iterations
for n=1:18
	res = fnlCg_xp(res,param);
	im_res = XFM'*res;
    if n==600
        ph=angle(im_res);
        param.FT=p2DFTxp(mask,N,cos(ph)+1i*sin(ph));
        im_res=im_res.*(cos(ph)-1i*sin(ph));
        res=XFM*im_res;
    end
    
    if ~nodisplay
        figure(101);
	 for i=1:length(isl)
      subplot(length(isl),2,(i-1)*2+1);imagesc(abs(im_res(:,:,isl(i))));colormap(gray);
      subplot(length(isl),2,(i-1)*2+2);imagesc(angle(im_res(:,:,isl(i))));
      drawnow;
     end
    end
end



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
param.TVWeight =TVWeight;%TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;%xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;


in=sum(data,4);
nm=sum(mask,4);
nm(nm==0)=1;
in=in./nm;
%im_dc=FT'*(repmat(in,[1,1,1,size(data,4)]));
im_mean=ifft2c(repmat(in,[1,1,1,size(data,4)]));
param.data = data-fft2c(im_mean);


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
%param.TV=Wavelet_ISU(sz(1:2));
param.TV=TVOPxp;
param.XFM=XFM;

if ~nodisplay
  figure(100);
  for i=1:length(isl)
    subplot(1,length(isl),i); imshow(abs(im_mean(:,:,isl(i))),[]); 
    drawnow;   
  end
      
 end

%res = param.XFM*im_dc;
res=param.XFM*zeros(size(im_mean));
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
    
    im_res=im_res+im_mean;
    if ~nodisplay
	  figure(100);
      for i=1:length(isl)
         subplot(1,length(isl),i); imshow(abs(im_res(:,:,isl(i))),[]);
         drawnow;   
      end
      %subplot(1,2,2);imshow(angle(im_res(:,:,isl)),[]); drawnow;
    end
end


function [sgn,pm]=reorient_varian2afni(orient)
        
        sgn=zeros(1,3);
        pm=zeros(1,3);
 
   
            
    switch orient
        case 'trans'
            order = {'A2P','R2L','S2I'};
        case 'trans90'
            order = {'L2R','A2P','S2I'};
        case 'cor'
            order = {'S2I','R2L','P2A'};
        case 'cor90'
            order = {'L2R','S2I','P2A'};
        case 'sag'
            order = {'S2I','P2A','L2R'};
        case 'sag90'
            order = {'A2P','S2I','L2R'};
        otherwise
            error('unknown orient');
    end
   
    %              axial;  axial90;  coronal; coronal90; sag;  sag90   
    %- to +    pro: A2P;    L2R;      S2I;     L2R;      S2I;  A2P;
    %- to +    ppe: R2L;    A2P;      R2L;     S2I;      P2A;  S2I;
    %- to +    pss: S2I;    S2I;      P2A;     P2A;      L2R;  L2R;
    
       for i=1:3 
        switch order{i}
            case 'R2L'
             sgn(i)=1;
             pm(1)=i; 
            case 'L2R'
              sgn(i)=-1;
              pm(1)=i;
            case 'A2P'
              sgn(i)=1;
              pm(2)=i;
            case 'P2A'
               sgn(i)=-1;
               pm(2)=i;
            case 'I2S'
               sgn(i)=1;
               pm(3)=i;
            case 'S2I'
              sgn(i)=-1;
              pm(3)=i;
            otherwise
                error('unknow directin');
                
        end
       end
       
       
       
     function [sgn,pm]=reorient_sdt2afni(orient)
        
        sgn=zeros(1,3);
        pm=zeros(1,3);
 
   if ~iscell(orient)
            
    switch orient
        case 'trans'
            order = {'A2P','R2L','S2I'};
        case 'trans90'
            order = {'L2R','A2P','S2I'};
        case 'cor'
            order = {'S2I','R2L','P2A'};
        case 'cor90'
            order = {'L2R','S2I','P2A'};
        case 'sag'
            order = {'S2I','P2A','L2R'};
        case 'sag90'
            order = {'A2P','S2I','L2R'};
        otherwise
            error('unknown orient');
    end
   
   else
       order=orient;
   end
    %L2R,A2P,S2I
    
       for i=1:3 
        switch order{i}
            case 'L2R'
             sgn(i)=1;
             pm(1)=i; 
            case 'R2L'
              sgn(i)=-1;
              pm(1)=i;
            case 'A2P'
              sgn(i)=1;
              pm(2)=i;
            case 'P2A'
               sgn(i)=-1;
               pm(2)=i;
            case 'S2I'
               sgn(i)=1;
               pm(3)=i;
            case 'I2S'
              sgn(i)=-1;
              pm(3)=i;
            otherwise
                error('unknow directin');
                
        end
       end                   
