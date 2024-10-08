function gems_csrecon(fid_prefix,method,ref,ref_ind,recon_ind)
%gems_csrecon(fid_prefix,method[,ref,ref_ind,recon_ind])

tstart=tic;
np=readPar(fid_prefix,'np');
nv=readPar(fid_prefix,'nv');
ns=readPar(fid_prefix,'ns');
npecs=readPar(fid_prefix,'npecs');
m2=mask_cs1d(fid_prefix);

pro=readPar(fid_prefix,'pro');
ppe=readPar(fid_prefix,'ppe');
lro=readPar(fid_prefix,'lro');
lpe=readPar(fid_prefix,'lpe');
%figure;imagesc(m2);drawnow;
arraydim=readPar(fid_prefix,'arraydim');
ne=arraydim;

z=read_fid([fid_prefix,'.fid']);

z=reshape(z,[np/2,ns,npecs,ne]);

z2=zeros(np/2,nv,ns,ne);
                                                                                                                                                                                               
for i=1:ns
    for j=1:ne  
     z2(:,m2(:,i,j)>0,i,j)=z(:,i,:,j);
    end
end
z2=fov_shift_kspace(z2,[0,-ppe],[lro,lpe]);
m2=shiftdim(m2,-1);
m2=repmat(m2,[np/2,1,1,1]);
orient = readPar(fid_prefix,'orient');
orient=orient(2:end-1);
fname=sprintf('%s_%s',fid_prefix,method);
switch method
    case 'nesta'
     
        recon=zeros(size(z2));
        fz2=ifft1c(z2,1);
        for i=1:size(z2,3)  
       
          for k=1:size(z2,1)            
          fprintf('Slice %d, Point %d\n',i,k);
        
        %  recon(k,:,i,:,j)=ch_scl(j)*ktFOCUSS_new(z2(k,:,i,:,j),[],[],3,m2(k,:,i,:),0.5,lambda,60,2,1);  
          recon(k,:,i,:)=ktNESTA_UP(fz2(k,:,i,:),m2(k,:,i,:));
        
          end
         
        end
         
        fprintf('total time = %f', toc);
        
    case 'ktf'
      %  z2=z2.*m2;
        recon=zeros(size(z2));
        
        for i=1:size(z2,3)            
          fprintf('Slice %d\n',i);
          recon(:,:,i,:)=ktFOCUSS_new(z2(:,:,i,:),[],[],3,m2(:,:,i,:),0.5,0.01,60,2,1);  
        end
    case 'lowk'
        a=repmat(sum(m2,4),[1,1,1,size(m2,4)]);
        z3=z2;
        z3(a<size(m2,4))=0;
        recon=ifft2c(z3);
      
    case 'omp'
         opts = [];
    opts.slowMode = 1;
    opts.printEvery     = 25;
    zmean=sum(z2,4)./sum(m2,4);
    zmean(sum(m2,4)==0)=0;
    z2res=z2-repmat(zmean,[1,1,1,size(m2,4)]);
     z2res=z2res(:,:,1,1:end/2);
     m2=m2(:,:,1,1:end/2);
     
        A=@(b) omp_mat(b,m2);
        At=@(b) omp_matt(b,m2);
        recon=OMP({A,At},z2res(m2>0),64*32,[],opts);
    disp('');         
        
  %   recon=flipdim(recon,2);
    case 'ktfklt'
       
        recon=zeros(size(z2));
        
        if ~exist([fid_prefix,'_ktf+orig.HEAD'],'file') || ~exist([fid_prefix,'_ktf_ph+orig.HEAD'],'file') 
          prev=zeros(size(z2));
          for i=1:size(z2,3)            
            fprintf('Slice %d\n',i);
            prev(:,:,i,:)=ktFOCUSS_new(z2(:,:,i,:),[],[],3,m2(:,:,i,:),0.5,0.01,60,2,1);  
   
          end
        else
            mag=BrikLoad([fid_prefix,'_ktf+orig.HEAD']);
            ph=BrikLoad([fid_prefix,'_ktf_ph+orig.HEAD']);
            prev=mag.*exp(1i*ph);
        end
            
        
        for i=1:size(z2,3)            
          fprintf('Slice %d\n',i);
          for iklt=1:3
           recon(:,:,i,:)=ktFOCUSS_KLT_new(z2(:,:,i,:),prev,[],[],m2(:,:,i,:),0.5,0.01,60,2,1);  
           prev=recon(:,:,i,:);
          end
        end
        
        %%
    case 'ktfxp'   
        z2=z2.*m2;
        recon=zeros(size(z2));
        for i=4%1:size(z2,3)
          fprintf('Slice %d\n',i);
          recon(:,:,i,:)=ktFOCUSS_xp(z2(:,:,i,:),m2(:,:,i,:),0.1,1/2);
        end
     recon=flipdim(recon,2);
    
    case 'lowrank'
        
        ref=load('refg_8_12_40_TR1.0_5cyc');
        ref2=[ones(1,length(ref));ref(:)'];
      xfm=tTransform(ref2);
      TVWeight=0;
      xfmWeight=0;
      z_init=sum(z2,4);
      sm=sum(m2,4);
      sm(sm==0)=1;
      z_init=z_init./sm;
      z_init=repmat(z_init,[1,1,1,size(m2,4)]);
      for i=1:size(z2,3)
       recon(:,:,i,:)=run_cs(z2(:,:,i,:),z_init(:,:,i,:),[],[],m2(:,:,i,:),xfm,TVWeight,xfmWeight,4,true);
      end
     case 'modcs'
         arraydim=readPar(ref,'arraydim');
         ne=arraydim;
         zref=read_fid([ref,'.fid']);
         zref=reshape(zref,[np/2,ns,nv,ne]);
         zref=mean(zref(:,:,:,ref_ind),4);
         zref=permute(zref,[1,3,2]);
         zref=fov_shift_kspace(zref,[-pro,-ppe],[lro,lpe]);
         z2=cat(4,zref,z2(:,:,:,recon_ind));
         mref=ones(size(zref));
         m2=cat(4,mref,m2(:,:,:,recon_ind));
         XFM=Wavelet_ISU([size(z2,1),size(z2,2)]);
         recon = modcsresCausalDetection_xp(z2,m2,XFM,1);
         recon=flipdim(recon,2);
         fid_prefix=[fid_prefix,'_3'];
     case 'single'
           
      for nt=1%1:size(z2,4)
          
          
         XFM=Wavelet_ISU([size(z2,1),size(z2,2)]);
         mu=10;
         lambda=mu;
         gamma=mu/1000;
         z2=z2/max(abs(z2(:)));
         
         recon(:,:,:,nt) = split_bregman(m2(:,:,:,nt),z2(:,:,:,nt), mu,[lambda,lambda*0.1],gamma, 1, 6,XFM,0.01,1);
         ph=angle(recon(:,:,1,1));
    %     param.FT=p2DFTxp(mask,N,cos(ph)+1i*sin(ph));
     %mu=100;
        % lambda=mu*0.01;
        % gamma=mu/1000;
        % split_bregman(m2(:,:,1,nt),z2(:,:,1,nt), mu, lambda, 0,gamma, 30,
        % 6,XFM,0.01,cos(ph)+1i*sin(ph));  %somehow making phase correction
        % lead to divergent behavior in split_bregman.
        
       %  run_cs_splitb(squeeze(z2(:,:,islice,nt)),m2,false);
      %      recon(:,:,:,nt)=run_cs(squeeze(z2(:,:,1,nt)),m2(:,:,1,nt),false,1);
           
      end
    case 'ft'
        XFMw=0.1;
        TVw=0.01;
        niter=10;
        xfm=FTt;
        z2_init=mean_sampled_kt(z2,m2);
        z2_init=repmat(z2_init,[1,1,1,size(m2,4)]);
        for i=1:size(z2,3)
          %recon(:,:,i,:)=run_cs_FTt(z2(:,:,i,:),m2(:,:,i,:),true,1,'ft',niter,XFMw,TVw);
          recon(:,:,i,:)=run_cs(z2(:,:,i,:),z2_init(:,:,i,:),[],[],m2(:,:,i,:),xfm,TVw,XFMw,false);
          
        end
        
    fname=sprintf('%s_%s_Wxfm%4.3f_TVw%4.3f_iter%d',fid_prefix,method,XFMw,TVw,niter);
    
    case 'ftdtr'
        
        XFMw=0.1;
        TVw=0.01;
        niter=10;
        xfm=FTt;
        z2_init=mean_sampled_kt(z2,m2);
        z2_init=repmat(z2_init,[1,1,1,size(m2,4)]);
        for i=1:size(z2,3)
          recon(:,:,i,:)=run_cs(z2(:,:,i,:)-z2_init(:,:,i,:),zeros(size(z2(:,:,1,:))),[],[],m2(:,:,i,:),xfm,TVw,XFMw,niter,false);
          
        end
        
    fname=sprintf('%s_%s_Wxfm%4.3f_TVw%4.3f_iter%d',fid_prefix,method,XFMw,TVw,niter);
    
    case 'klt'
        recon=run_cs_FTt(z2,m2,false,1:9,'klt');
        

      otherwise 
          error('unknow methods\n');
end

    
switch orient 
   case 'trans90'  %1,2,3 in b - l2r,u2d in stimulate - l2r, a2p,s2i in vnmrj
   
   recon=flipdim(recon,2);  
%   b=flipdim(b,1);  
  order = {'L2R','A2P','S2I'};
   
  case 'trans'
   recon = permute(recon,[2,1,3,4]);   %do I really need this? yes 
   recon=flipdim(recon,2); 
   order={'L2R','A2P','S2I'};
   
 case 'sag90'  % 1,2 in b - l2r, u2d in stimulate - a2p, s2i in vnmrj
  %  b = permute(b,[2,1,3,4]);
      recon=flipdim(recon,1);    
      recon=flipdim(recon,2);
      order = {'A2P','S2I','L2R'};
      
    case 'sag' % 1,2 in b - l2r, u2d in stimulate - a2p, s2i in vnmrj
      recon=permute(recon,[2,1,3,4]);   
      recon=flipdim(recon,2); 
      order = {'A2P','S2I','L2R'};
      
    case 'cor90'  
      recon=  flipdim(recon,1);
      recon=flipdim(recon,2);
      order = {'L2R','S2I','P2A'};
    case 'cor'
      recon=permute(recon,[2,1,3,4]);
      recon=flipdim(recon,2);
      
      order = {'L2R','S2I','P2A'};
      
    otherwise
        
end
%{
tmp1=rdSdt([fid_prefix,'_',method]);
tmp2=rdSdt([fid_prefix,'_',method,'_ph']);
recon=tmp1.*exp(1i*tmp2);
%}

%{
writesdt4(abs(recon),[fid_prefix,'_',method]);
writesdt4(angle(recon),[fid_prefix,'_',method,'_ph']);
%}

[sgn,pm]=reorient_sdt2afni(order);

tmp=recon;
for i=1:3
    if sgn(i)==-1
        tmp=flipdim(tmp,i);
    end
end

afni=permute(tmp,[pm,4]);
    

%%
d=fid_prefix;
lro = readPar(d,'lro');%in cm
lpe = readPar(d,'lpe');%in cm
thk = readPar(d,'thk');%in mm
tr = readPar(d,'tr');
pss=readPar(d,'pss'); %in cm
%img = slice_reorder(img,pss);    

    pro=readPar(d,'pro'); %in cm
    ppe=readPar(d,'ppe'); %in cm
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
    np=readPar(fid_prefix,'np');
    nread=np/2;
    nphase=readPar(fid_prefix,'nv');
    delta = [lro*10/nread,lpe*10/nphase,thk];
  
    center= center.*sgn;
    center=center(pm);
    
    delta=delta(pm);
    delta([1,3])=-delta([1,3]);
    
    if length(sz)==2
        sz(3)=1;
    end
    %    
  %   [img,orig_delta] = reorient_data(img,orig_delta,orient(2:end-1));
  %   should not need to reorient if data is from sdt file.  6/4/2012
    info.ORIGIN = center-delta.*(sz(1:3)/2-[0.5,0.5,0.5]);
    info.DELTA = delta;
    info.TAXIS_FLOATS = [0,tr,0,0,0];
    info.ORIENT_SPECIFIC = [1 3 5]; %L2R,A2P,S2I  % afni convention -2+ is R2L, A2P, I2S
    write_afni(abs(afni),fname,info);
    write_afni(angle(afni),[fname,'_ph'],info);

    fprintf('gems_csrecon finished in %4.1f s\n',toc(tstart));
    

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


    
function [sgn,pm]=reorient_varian2afni(orient)
    % coordinate conversion.  
    % covert it to afni convention of -2+ is R2L, A2P, I2S
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
    % data conversion
    % orient is the direction of the original data
    % pm is the matrix used to convert it to afni order of [1,3,5] or L2R,A2P,S2I
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
                

    function out=omp_mat(b,m)
        
        out=fft2c(ifft1c(b,4)).*m; 
        out=out(:);
        
        
   function out=omp_matt(b,m)
            
        d=zeros(size(m));
        d(m>0)=b;
        out=ifft2c(fft1c(d,4));
        out=out(:);        
        
        

