function gems_csreconbr(a,format)
%gems_csrecon(a,reduction,format)

if ~exist('format','var')
    format='a';
end
tstart=tic;

of1=readbPar(fullfile(a,'method'),'PVM_SPackArrPhase1Offset');
of0=readbPar(fullfile(a,'method'),'PVM_SPackArrReadOffset');
fov=readbPar(fullfile(a,'method'),'PVM_Fov');

ns=readbPar(fullfile(a,'method'),'PVM_SPackArrNSlices');
obj=readbPar(fullfile(a,'method'),'PVM_ObjOrderList');


%a=fov_shift_kspace(a3,[0,-of1],fov);
reduction = readbPar(fullfile(a,'method'),'R');
d=readb_fid(a);

sw=readbPar(fullfile(a,'method'),'PVM_EffSWh');
D=readbPar(fullfile(a,'acqp'),'D');
pe=readbPar(fullfile(a,'acqp'),'ACQ_spatial_phase_1');
NR = readbPar(fullfile(a,'method'),'PVM_NRepetitions');
kzero = find(pe==0);
riset=D(4);

navgrad=readbPar(fullfile(a,'method'),'NavGradOn',false);

rampt = D(5);%readbPar(fullfile(a,'method'),'PVM_RampTime');



if strcmp(navgrad,'No')
    askip=0;
else
    askip = round((riset+rampt)*sw);
end

mtrx=readbPar(fullfile(a,'method'),'PVM_Matrix');
pe=round(pe*mtrx(2)/2)+mtrx(2)/2+1;
nav=d(askip+1:askip+mtrx(1),:,:,:);
data=d(end-mtrx(1)+1:end,:,:,:);
data2=data;
[nch,ch_scl]=nChan(a);

if mod(NR,length(pe)/mtrx(2))~=0
   warning('mod(NR,length(pe)/mtrx(2))~=0');
end




iNR = length(pe)/mtrx(2);
%data=reshape(data,[mtrx(1),nch,ns,mtrx(2)/reduction,reduction*iNR,NR/iNR]);
data=reshape(data,[mtrx(1),nch,ns,mtrx(2)/reduction,reduction*NR]);
nav=reshape(nav,[mtrx(1),nch,ns,mtrx(2)/reduction,reduction*NR]);

pe=reshape(pe,[mtrx(2)/reduction,reduction*iNR]);

data=permute(data,[1,4,3,5,2]);
nav=permute(nav,[1,4,3,5,2]);

%data=phaseCorr(data,nav(1,:,:,:,:));

figure;
nav1=nav(:,:,1,:,1);
data1=data(:,:,1,:,1);

plot(abs(nav1(:,kzero(1))));
hold on;
y=flipdim(data1(:,kzero(1)),1);
plot(abs(y),'r');
legend('Nav Echo','k=0 data');

z2=zeros([mtrx',ns,reduction*NR,nch]);
m2=zeros(mtrx(1),mtrx(2),ns,reduction*NR);

for j=1:NR*reduction
  jj=mod(j-1,iNR*reduction)+1;
  z2(:,pe(:,jj),:,j,:)=data(:,:,:,j,:);
  m2(:,pe(:,jj),:,j)=1;
end


%%
of1=readbPar(fullfile(a,'method'),'PVM_SPackArrPhase1Offset');
of0=readbPar(fullfile(a,'method'),'PVM_SPackArrReadOffset');
fov=readbPar(fullfile(a,'method'),'PVM_Fov');
z2=fov_shift_kspace(z2,[0,-of1],fov);



%%
lambda=0;
method='nesta';
switch method
    case 'nesta'
        tic;
        recon=zeros(size(z2));
        fz2=ifft1c(z2,1);
        for i=1:size(z2,3)  
       
        for j=1:nch
          for k=1:size(z2,1)            
          fprintf('Channel %d, Slice %d, Point %d\n',j,i,k);
        
        %  recon(k,:,i,:,j)=ch_scl(j)*ktFOCUSS_new(z2(k,:,i,:,j),[],[],3,m2(k,:,i,:),0.5,lambda,60,2,1);  
          recon(k,:,i,:,j)=ktNESTA_UP(fz2(k,:,i,:,j),m2(k,:,i,:));
        
          end
         write_afni(abs(recon(:,:,i,:,j)),sprintf('scan8_slice%d_ch%d',i,j));
        
       end
        
        end
         
        fprintf('total time = %f', toc);
        
      %  save gems_csreconbr
    case 'ktf'
      %  z2=z2.*m2;
        recon=zeros(size(z2));
        
        
        for j=1:nch
         for i=1:size(z2,3)            
          fprintf('Channel %d, Slice %d\n',j,i);
          recon(:,:,i,:,j)=ch_scl(j)*ktFOCUSS_new(z2(:,:,i,:,j),[],[],3,m2(:,:,i,:),0.5,lambda,60,2,1);  
         end
        end
    case 'lowk'
        a=repmat(sum(m2,4),[1,1,1,size(m2,4)]);
        z3=z2;
        z3(a<size(m2,4))=0;
        recon=ifft2c(z3);
    case 'omp'
        
        
        A=@(b) omp_mat(b,m2);
        At=@(b) omp_matt(b,m2);
         opts = [];
    opts.slowMode = 1;
    opts.printEvery     = 25;
    zmean=sum(z2,4)./sum(m2,4);
    zmean(sum(m2,4)==0)=0;
    z2res=z2-repmat(zmean,[1,1,1,size(m2,4)]);
        
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
           recon(:,:,i,:)=ktFOCUSS_KLT_new(z2(:,:,i,:),prev,[],[],m2(:,:,i,:),0.5,0.01,100,2,1);  
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

sz=size(recon);
recon=reshape(recon,[sz(1:3),sz(4)*sz(5)]);
   

recon(:,:,obj+1,:)=recon;


if format=='a'
 write_afni(abs(recon),[a,'_',method,'_nonav_',num2str(lambda)]);
elseif format=='s'
    writesdt4(abs(recon),[a,'_',method]);
else
    error('unknown format');
end

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

       