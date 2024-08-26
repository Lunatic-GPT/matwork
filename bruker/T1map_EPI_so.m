function T1map_EPI_so(scand,fM0,mask)

tic;
if ~exist(sprintf('%s_recon+orig.HEAD',scand),'file');
epiShaperbr(scand,true);
end
if ~exist(sprintf('%s_recon+orig.HEAD',fM0),'file');
epiShaperbr(fM0,true);
end

a=BrikLoad(sprintf('%s_reconchan+orig.HEAD',scand));
ph=BrikLoad(sprintf('%s_reconph+orig',scand));

M0=BrikLoad(sprintf('%s_reconchan+orig.HEAD',fM0));
M0ph=BrikLoad(sprintf('%s_reconph+orig',fM0));

M0=(M0.*exp(1i*M0ph));

ph=reshape(ph,[size(M0,1),size(M0,2),size(M0,3),size(ph,4)/4,4]);

M0=reshape(M0,[size(M0,1),size(M0,2),size(M0,3),size(M0,4)/4,4]);
M0=(mean(M0,4));
M0ph=repmat(angle(M0),[1,1,1,size(ph,4),1]);

sgn=abs(M0ph-ph)<pi/2 | abs(M0ph-ph)>pi*3/2;

sgn= (2*(sum(sgn,5)>2))-1;

a=reshape(a,size(ph));

a=sqrt(mean(abs(a).^2,5));
M0=sqrt(mean(abs(M0).^2,5));

M0=repmat(M0,[1,1,1,size(a,4)]);
a=1-a.*sgn./M0;

if exist('mask','var')
mask=load(mask);
mask=mask.roi;
else
    mask=ones(size(a,1),size(a,2),size(a,3));
end
ideg=[49,65,1];

    debug=false;
%M0=
TRimg=readbPar([scand,'/method'],'TRimg');

ns=readbPar([scand,'/acqp'],'NSLICES');

freq=readbPar([scand,'/acqp'],'ACQ_O1_list');

firstTI=readbPar([scand,'/method'],'TI0');

nTR=readbPar([scand,'/method'],'PVM_NRepetitions');

TI=zeros(nTR,ns);


freqlist=sort(unique(freq),'ascend');

%% assuming the first slice has the lowest frequency.

for i=1:nTR
    for j=1:ns
     ind=mod(j+ns*(i-1)-1,length(freq))+1;
     
     ind2=find(freq(ind)==freqlist);
     TI(i,ind2)=firstTI+TRimg/ns*(j-1);
    end
end
options=optimset('MaxIter',30,'Display','off','FunValCheck','off','TolFun',0.001,'TolX',0.001);
t1=zeros(size(a,1),size(a,2),size(a,3));
res=t1;
if ~debug
matlabpool(4);
end

parfor i=1:size(a,1)
    disp(i);
    t2=zeros(size(a,2),size(a,3));
    res2=t2;
    for j=1:size(a,2)
        for k=1:size(a,3)
           if mask(i,j,k)==0
             continue;
           end
            y = squeeze(a(i,j,k,:));
            te2=TI(2:end,k)/1000;
               
            if ~debug
              [b,r]=lsqcurvefit(@exp_decay,[max(y),2],te2(:),y(:),[],[],options);
               
            if ~isnan(b(1)) && ~isnan(b(2))
             ss = sum((y-mean(y)).^2);   
             t2(j,k) = 1/b(2);
             res2(j,k) = 1-sum(r.^2)/ss;
          %   beta(i,j,k,ir)=b(3);
            end
            
            else
                
            if i==ideg(1)&&j==ideg(2)&&k==ideg(3)
             [b,r]=lsqcurvefit(@exp_decay,[max(y),2],te2(:),y(:),[],[],options);
             disp(b);
              figure;plot(te2,y,'o');
              hold on;
              plot(te2,exp_decay(b,te2),'r-');
                
            end
            
            end         
        end
    end
   res(i,:,:)=res2;
   t1(i,:,:)=t2;
end

if ~debug
matlabpool close;
end
if ~debug
write_afni(t1,sprintf('%s_T1map_expdecay',scand));
end

fprintf('T1map_EPI_so finished in %f s\n',toc);


function epiShaperbr(fid_dir,doref)
% epiShaperbr(fid_dir,doref)

if ~exist('doref','var')
    doref=true;
end


a=readb_fid(fid_dir);

nro=size(a,1);
sz=readbPar(fullfile(fid_dir,'method'),'PVM_EncMatrix');

mt=readbPar(fullfile(fid_dir,'method'),'PVM_Matrix');

ns=readbPar(fullfile(fid_dir,'acqp'),'NSLICES');

nseg=readbPar(fullfile(fid_dir,'method'),'NSegments');
NR=readbPar(fullfile(fid_dir,'acqp'),'NR');

acqmod=readbPar(fullfile(fid_dir,'acqp'),'ACQ_experiment_mode',false);
necho= readbPar(fullfile(fid_dir,'method'),'PVM_NEchoImages');

if ~strcmp('SingleExperiment',acqmod)
    rcvrs=readbPar(fullfile(fid_dir,'acqp'),'ACQ_ReceiverSelect',false);
    nyes=strmatch('Yes',rcvrs);
    nch=length(nyes);    
   
else
    nch=1;
end


a=reshape(a,[sz(1),sz(2)/nseg,nch,necho,ns,nseg,NR]);

a=permute(a,[1,2,6,5,4,7,3]);  %sz(1),sz(2)/nseg,nseg,ns,necho,NR,nch

a=reshape(a,[sz(1),sz(2),ns,necho,NR,nch]);
a = convertTraces(a);  % reverse the direction of even number traces.
a2=zeros(size(a));

freq=readbPar([fid_dir,'/acqp'],'ACQ_O1_list');
firstTI=readbPar([fid_dir,'/method'],'TI0');
nTR=readbPar([fid_dir,'/method'],'PVM_NRepetitions');
freqlist=sort(unique(freq),'ascend');

for i=1:nTR
    for j=1:ns
     ind=mod(j+ns*(i-1)-1,length(freq))+1;
     
     ind2=find(freq(ind)==freqlist);
     a2(:,:,ind2,:,i,:)=a(:,:,j,:,i,:);
    end
end

a=a2;
figure;
for i=1:nch
   subplot(nch,1,i); 
    imshow(abs(a(:,:,1,1,1,i)),[]);
end

ky=readbPar(fullfile(fid_dir,'method'),'PVM_EncSteps1');

a2=a;
a(:,ky-min(ky)+1,:,:,:,:)=a2;
%a=reshape(a,[sz(1),sz(2),ns,length(a(:))/sz(1)/sz(2)/ns/nch/necho,nch]);

ref=a(:,:,:,:,1,:);

tmp=reshape(ref,[size(a,1)*size(a,2)*size(a,3)*size(a,4),size(a,6)]);
[tmp,ind]=max(max(abs(tmp),[],1),[],2);

if doref
    a = fidPhaseCorr2(a(:,:,:,:,2:end,:),ref(:,:,:,1,1,ind));
else
   % a=a(:,:,:,:,2:end,:);
end
sz=size(a);
a=reshape(a,[sz(1:3),sz(4)*sz(5),size(a,6)]);
if mt(1)<sz(1)
    
    kx=get_kxtraj(fid_dir);
    
%kx=readbPar(fullfile(fid_dir,'method'),'PVM_EpiTrajAdjkx');
%{
kx=[0 -0.0818895161648283 -0.117473651411874 -0.0848028534241132 
0.0339713879037591 0.233754350190432 0.517501461981041 0.887466000227911 
1.34196690472017 1.87092704371424 2.48274936276364 3.15944118378767 
3.89762551480635 4.68551478273991 5.52330660620812 6.38716747354377 
7.2909643964452 8.19913518935879 9.12890992180585 10.0682239216206 
11.0237367268672 11.9867021304038 12.9449003331803 13.912639003131 
14.8650856769237 15.8372357116181 16.8214978754329 17.7951326521657 
18.7769835508448 19.7544611659737 20.748701406813 21.7409863990848 
22.7349672194854 23.7281132320737 24.7205331490342 25.7256820771653 
26.7217065273873 27.7218010699124 28.7272127199198 29.7255625878683 
30.7328191323387 31.7239637436798 32.7192874765882 33.7156352237497 
34.7160074957933 35.7091396982144 36.7037086114631 37.6960502693826 
38.6673414541929 39.6445001118114 40.598184110452 41.5282510455021 
42.4314957990014 43.3039280740382 44.1222317850021 44.8958901048685 
45.6126336197549 46.2659706958736 46.8575997908682 47.3822753945902 
47.8271268317848 48.1928583840579 48.4877209197763 48.7075611696063];
%b0=readbPar(fullfile(fid_dir,'method'),'PVM_EpiTrajAdjb0');
kx=kx';
kx=kx(:)';
%}
while 1
    if kx(end)==0
        kx(end)=[];
 %       b0(end)=[];
    else
        break;
    end
end


fprintf('Ramp sampling correction on.\n');

if length(kx)<sz(1)
kx2=interp1((0:length(kx)-1)*(sz(1)-1)/(length(kx)-1),kx*(sz(1)-1)/(length(kx)-1),0:sz(1)-1);
%b02=interp1((0:length(kx)-1)*(sz(1)-1)/(mt(1)-1),b0*(sz(1)-1)/(mt(1)-1),0:sz(1)-1);
else
    kx2=kx;
end

sz=size(a);
sz(1)=mt(1);
a2=zeros(sz);
for i=1:size(a,2)
    for j=1:size(a,3)
        for k=1:size(a,4)
            for l=1:size(a,5)
      %     a2(:,i,j,k)=interp1(kx2,a(:,i,j,k),0:mt(1)-1);
           a2(:,i,j,k,l)=interp1(kx2,a(:,i,j,k,l),-mt(1)/2:mt(1)/2-1);
            end
        end
    end
end
else
    fprintf('Ramp sampling correction off.\n')
    a2=a;
end

a2(isnan(a2))=0;
%% shift image
of1=readbPar(fullfile(fid_dir,'method'),'PVM_SPackArrPhase1Offset');
of0=readbPar(fullfile(fid_dir,'method'),'PVM_SPackArrReadOffset');
fov=readbPar(fullfile(fid_dir,'method'),'PVM_Fov');

%z2=fov_shift_kspace(a2,[-of0,-of1],fov);

%z2=a;
if mt(2)>size(z2,2)
  nneg=readbPar(fullfile(fid_dir,'method'),'PVM_EncCentralStep1');
  img=zeros([mt',size(z2,3),size(z2,4),size(z2,5)]);
  for i=1:size(z2,5)
  img(:,:,:,:,i)=partialFT(z2(:,:,:,:,i),nneg-1);
  
  end
  ph=angle(img);
  ph=reshape(ph,[size(ph,1),size(ph,2),size(ph,3),size(ph,4)*size(ph,5)]);
  img2=sqrt(mean(abs(img).^2,5));
  img=reshape(img,[size(ph,1),size(ph,2),size(ph,3),size(ph,4)*size(ph,5)]);
else
 fz=fft2c(z2);
 fz=flipdim(fz,1);
 fz=flipdim(fz,2);
 
  ph=angle(fz);
  ph=reshape(ph,[size(ph,1),size(ph,2),size(ph,3),size(ph,4)*size(ph,5)]);
  img=reshape(fz,[size(ph,1),size(ph,2),size(ph,3),size(ph,4)*size(ph,5)]);
  
  img2=sqrt(mean(abs(fz).^2,5));
end

%objlist=readbPar([fid_dir,'/method'],'PVM_ObjOrderList');
%img(:,:,objlist+1,:)=img;
if doref
 write_afni(abs(img),[fid_dir,'_reconchan']);
 write_afni(ph,[fid_dir,'_reconph']);
 write_afni(img2,[fid_dir,'_recon']);
 
 %write_afni(angle(fz),[fid_dir,'_reconph']);
else   
 write_afni(img2,[fid_dir,'_noref']);
 write_afni(abs(img),[fid_dir,'_norefchan']);
 write_afni(ph,[fid_dir,'_norefph']);
 %write_afni(angle(fz),[fid_dir,'_ph_noref']);
end
    
function z = convertTraces(z_tmp)
 % reverse the direction of even number traces
 z = zeros(size(z_tmp));
 
 z(:,1:2:end-1,:,:,:,:,:) = z_tmp(:,1:2:end-1,:,:,:,:,:);
 z(1:end,2:2:end,:,:,:,:,:) = z_tmp(end:-1:1,2:2:end,:,:,:,:,:);
        
        
   function z=fidPhaseCorr2(z,zref)
  %
        sz = size(z);
        fz = fft1c(z,1);
        fzref = fft1c(zref,1);
        sz_rep = sz;
        if length(sz)>=4 
          sz_rep(1:end-length(sz)+3)=1;
          fzref = repmat(fzref,sz_rep);
        end
        
        H=mean(abs(fz(:)));
        ph=conj(fz).*fzref./abs(fz)./abs(fzref);
     %   ph(:,1:2:end-1,:,:)=1;
        %ph(abs(fz)<0.1*H & abs(fzref)<0.1*H)=1;
        
          fzref = fzref./abs(fzref); 
          
        fz = fz.*conj(fzref);
        
        z=ifft1c(fz,1);
        
function kx=get_kxtraj(fid)
           
    fid=fullfile(fid,'method');
g=readbPar(fid,'PVM_EpiReadOddGrad');%0.326530612244898

g2=readbPar(fid,'PVM_EpiReadEvenGrad');%0.326530612244898

pt=readbPar(fid,'PVM_EpiPlateau');%0.1548
rt=readbPar(fid,'PVM_EpiRampTime'); %0.05
es=readbPar(fid,'PVM_EpiEchoSpacing');  %0.2548
const=readbPar(fid,'PVM_GradCalConst');  %31250 Hz/mm 
enc=readbPar(fid,'PVM_EncMatrix'); %91*64
sw=readbPar(fid,'PVM_DigSw');  %sw=g*fov*const

fov=readbPar(fid,'PVM_Fov');  %35*35

t=linspace(0,es,enc(1)+1);


y=zeros(1,enc(1));

for i=1:enc(1)
    if t(i)<rt
        y(i)=t(i)/rt;
    elseif t(i)<rt+pt
        y(i)=1;
    else
        y(i)=(es-t(i))/rt;
    end
end


kx=cumsum(y);
kx=kx-sum(y(1:round((enc(1)+1)/2)));
