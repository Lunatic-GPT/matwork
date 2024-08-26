function kdata_epi(fid_dir,doref)
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

isnav=readbPar(fullfile(fid_dir,'method'),'PVM_EpiPrefixNavYes',false);

acqmod=readbPar(fullfile(fid_dir,'acqp'),'ACQ_experiment_mode',false);
if ~strcmp('SingleExperiment',acqmod)
    rcvrs=readbPar(fullfile(fid_dir,'acqp'),'ACQ_ReceiverSelect',false);
    nyes=strmatch('Yes',rcvrs);
    nch=length(nyes);    
   
else
    nch=1;
end

if strcmp(isnav,'Yes')
  n_nav=readbPar(fullfile(fid_dir,'method'),'PVM_EpiPrefixNavSize');
  npt=sz(1)*sz(2)/nseg; 
  nav=a(1:n_nav,:,:);
  a=a(end-npt+1:end,:,:);
  a=reshape(a,[sz(1),sz(2)/nseg,ns,nseg,length(a(:))/sz(1)/sz(2)/ns]);

  nav=reshape(nav,[n_nav/2,2,ns,nseg,length(nav(:))/n_nav/ns/nseg]);
  
  ind=1;
  for i=1:size(nav,3)
      for j=1:size(nav,4)
          for k=1:size(nav,5)
              x=nav(:,:,i,j,k);
               p = polyfit((1:length(x(:)))',angle(x(:)),1);
               y = a(:,:,i,j,k);
               xx=(nro-length(y(:))+1):nro;
               y=y(:).*exp(-1i*(p(1)*xx'+p(2)));
             % pp=[-2.5113,-2.501,-2.5271,-2.4885,-2.535,-2.5166];
          %    pp=[-2.5,-2.5,-2.5,-2.5,-2.5,-2.5];
           %  pp=zeros(1,6);
              %y=transpose(y(:)).*exp(-1i*(pp(ind)));       
             % y=y(:);
              a(:,:,i,j,k)=reshape(y,[size(a,1),size(a,2)]);
        ind=ind+1;       
              disp(p);
          end
      end
  end
  %{
  ang=repmat(nav(:,1,:,:,:),[1,size(a,2)/2,1,1,1]);
  ang=angle(ang);
  a(:,1:2:end-1,:,:,:)=a(:,1:2:end-1,:,:,:).*exp(-1i*ang);
  
  
  ang=repmat(nav(:,2,:,:,:),[1,size(a,2)/2,1,1,1]);
  ang=angle(ang);
  a(:,2:2:end,:,:,:)=a(:,2:2:end,:,:,:).*exp(-1i*ang);
  %}
end

a=reshape(a,[sz(1),sz(2)/nseg,nch,ns,nseg,length(a(:))/sz(1)/sz(2)/ns/nch]);

a=permute(a,[1,2,5,4,6,3]);
a=reshape(a,[sz(1),sz(2),ns,length(a(:))/sz(1)/sz(2)/ns/nch,nch]);
a = convertTraces(a);  % reverse the direction of even number traces.

%%
% appears to be no improvement of image quality with b0 correction. 
%{  
b0=readbPar(fullfile(fid_dir,'method'),'PVM_EpiTrajAdjb0');

while 1
    if b0(end)==0
        b0(end)=[];
    else
        break;
    end
end

b02=interp1((0:length(b0)-1)*(sz(1)-1)/(mt(1)-1),b0*(sz(1)-1)/(mt(1)-1),0:sz(1)-1);

b02=repmat(b02',[1,size(a,2),size(a,3),size(a,4)]);
a=a.*exp(-1i*b02);
%}

%%
figure;
for i=1:nch
   subplot(nch,1,i); 
    imshow(abs(a(:,:,1,1,i)),[]);
end

ky=readbPar(fullfile(fid_dir,'method'),'PVM_EncSteps1');
a2=a;
a(:,ky-min(ky)+1,:,:,:)=a2;

ref=a(:,:,:,1,:);

tmp=reshape(ref,[size(a,1)*size(a,2)*size(a,3),size(a,5)]);
[tmp,ind]=max(max(abs(tmp),[],1),[],2);

if doref
    a = fidPhaseCorr2(a(:,:,:,2:end,:),ref(:,:,:,1,ind));
else
    a=a(:,:,:,2:end,:);
end

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

z2=fov_shift_kspace(a2,[-of0,-of1],fov);

save(sprintf('kdata_%s',fid_dir),'z2');

    
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
kx=kx-sum(y(1:(enc(1)+1)/2));




        
        
        
        




           
           
           
        
        