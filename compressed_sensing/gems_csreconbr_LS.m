function gems_csreconbr_LS(a,reduction,format)
%gems_csrecon(a,reduction,format)

if ~exist('format','var')
    format='s';
end
tstart=tic;


of1=readbPar(fullfile(a,'method'),'PVM_SPackArrPhase1Offset');
of0=readbPar(fullfile(a,'method'),'PVM_SPackArrReadOffset');
fov=readbPar(fullfile(a,'method'),'PVM_Fov');


%a=fov_shift_kspace(a3,[0,-of1],fov);

d=readb_fid(a);

sw=readbPar(fullfile(a,'method'),'PVM_EffSWh');
D=readbPar(fullfile(a,'acqp'),'D');
pe=readbPar(fullfile(a,'acqp'),'ACQ_spatial_phase_1');
NR = readbPar(fullfile(a,'method'),'PVM_NRepetitions');
kzero = find(pe==0);
riset=D(4);

rampt = D(5);%readbPar(fullfile(a,'method'),'PVM_RampTime');

askip = round((riset+rampt)*sw);
mtrx=readbPar(fullfile(a,'method'),'PVM_Matrix');
pe=round(pe*mtrx(1)/2)+mtrx(1)/2+1;
nav=d(1,:,:,:);

data=d(end-mtrx(1)+1:end,:,:,:);
data=phaseCorr(data,nav);

%{
figure;
plot(abs(nav(:,kzero(1))));
hold on;
y=flipdim(data2(:,kzero(1)),1);

plot(abs(y),'r');

legend('Nav Echo','k=0 data');

%}
if mod(NR,length(pe)/mtrx(2))~=0
  % error('mod(NR,length(pe)/mtrx(2))~=0');
end

iNR = length(pe)/mtrx(2);
data=reshape(data,[mtrx(1),iNR*reduction,mtrx(2)/reduction]);

pe=reshape(pe,[iNR*reduction,mtrx(2)/reduction]);
pe=pe';
data=permute(data,[1,3,2,4]);


sz=size(data);

d2=zeros([mtrx',iNR*reduction]);
m2=zeros(mtrx(1),mtrx(2),1,iNR*reduction);

 for i=1:iNR*reduction
  d2(:,pe(:,i),i)=data(:,:,i);
  m2(:,pe(:,i),1,i)=1;
 end

z2=reshape(d2,[mtrx',1,iNR*reduction]);

m2=reshape(m2,[mtrx',1,iNR*reduction]);
%%
of1=readbPar(fullfile(a,'method'),'PVM_SPackArrPhase1Offset');
of0=readbPar(fullfile(a,'method'),'PVM_SPackArrReadOffset');
fov=readbPar(fullfile(a,'method'),'PVM_Fov');
z2=fov_shift_kspace(z2,[0,-of1],fov);



%%
method='ktf';
switch method
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

    
if format=='a'
 write_afni(abs(recon),[a,'_',method]);
elseif format=='s'
    writesdt4(abs(recon),[a,'_',method]);
else
    error('unknown format');
end

    fprintf('gems_csrecon finished in %4.1f s\n',toc(tstart));
    
        
function d2=phaseCorr(d2,nav)

  ph=angle(nav);
  
  ph=repmat(ph,[size(d2,1),1,1,1]);
  
 d2=d2.*exp(-1i*ph);


   
