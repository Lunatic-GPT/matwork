function recon_UTE3d_Bruker(fname)

a=readb_fid(fname);

%trj=readTraj(fname);
%trj=readTraj('12');
trj=readTraj(fname);

mtrx=readbPar([fname,'/method'],'PVM_Matrix');
npro=readbPar([fname,'/method'],'NPro');

% do the following if a DCF file is not available
calc_DCF = false; % calculate density conpensation? 
if calc_DCF
  [DCF,gdata2]=calc_DCF(trj,mtrx(2));
  tmp.DCF=DCF;
else
 tmp=load('4_DCF_N32.mat');
end

sz=size(squeeze(a));
data=zeros([2,sz]);

data(1,:,:,:,:)=real(a);
data(2,:,:,:,:)=imag(a);

N=mtrx(1)*1.5;
% 
%  ib=1;
%  gdata=grid3_MAT_xp(data(:,ib:end,:),trj(:,ib:end,:),tmp.DCF(ib:end,:),N);
%  gdata=gdata(1,:,:,:)+1i*gdata(2,:,:,:);

k_not = [0.0, 0.0, 0.0];
DCF   = [1.0];

[x,y,z]=meshgrid(linspace(-0.5,0.5,N),linspace(-0.5,0.5,N),linspace(-0.5,0.5,N));
g2=zeros(3,N,N,N);
g2(1,:,:,:)=x;
g2(2,:,:,:)=y;
g2(3,:,:,:)=z;

rokern = grid3_MAT_2grid_real(1.0,k_not',g2,N);
rokern=ifft2c(rokern);
rokern=ifft1c(rokern,3);

gdata1=grid3_MAT_2grid_real(tmp.DCF.*squeeze(data(1,:,:,:)),trj,g2,N);
gdata2=grid3_MAT_2grid_real(tmp.DCF.*squeeze(data(2,:,:,:)),trj,g2,N);
gdata=gdata1+1i*gdata2;
img=ifft2c(squeeze(gdata));
img=ifft1c(img,3);

%gdata1=grid3_MAT_2grid_real(tmp.DCF.*squeeze(data(1,:,:,:)),trj,g2,32);
img2=img;
rokern=rokern/max(abs(rokern(:)));
mask=abs(rokern)>0.02;
img2(mask)=img(mask)./rokern(mask);
img2(~mask)=0;
img_abs=abs(img2);

save([fname,'_recon'],'img_abs');
%write_afni(abs(img2),[fname,'_recon']);
%{
sb=round((N-mtrx(1))/2)+1;
se=sb+mtrx(1)-1;
img2=img2(sb:se,sb:se,sb:se);

figure;
for i=1:32
subplot(4,8,i);
imshow(abs(img2(:,:,i)),[0,max(abs(img2(:)))]);
end
%}
function a=readTraj(fname)

 fid=fopen([fname,'/traj'],'r','ieee-le');
b=fread(fid,'double');

mtrx=readbPar([fname,'/method'],'PVM_Matrix');
npro=readbPar([fname,'/method'],'NPro');

fclose(fid);

a=reshape(b,[3,mtrx(1)/2,npro]);


function [res,nzero]=readb_fid(d)

if ~exist('d','var')
    d=pwd;
end


sz=readbPar(fullfile(d,'acqp'),'ACQ_size',true);

nr=readbPar(fullfile(d,'acqp'),'NR',true);
ni=readbPar(fullfile(d,'acqp'),'NI',true);

if numel(sz)>1
 ns = readbPar(fullfile(d,'acqp'),'NSLICES',true);
else
    ns=1;
end

fid=fopen(fullfile(d,'fid'),'r','ieee-le');

fmt=readbPar(fullfile(d,'acqp'),'GO_raw_data_format',false);

if strcmp(fmt,'GO_32BIT_SGN_INT')
 res=fread(fid,'int32');
 
else
    fprintf('%s :',fmt);
    fclose(fid);
    error('unknown format');
end

acqmod=readbPar(fullfile(d,'acqp'),'ACQ_experiment_mode',false);
if ~strcmp('SingleExperiment',acqmod)
    rcvrs=readbPar(fullfile(d,'acqp'),'ACQ_ReceiverSelect',false);
    nyes=strmatch('Yes',rcvrs);
    nch=length(nyes);    
   
else
    nch=1;
end

bs=readbPar(fullfile(d,'acqp'),'GO_block_size',false);

if strcmp('Standard_KBlock_Format',bs);
    sz1_old=sz(1)*nch;
 
 %s1=2.^ceil(log2(sz(1)*nch));
 %s2=2.^ceil(log2(sz(1)*nch/10))*10;
 %sz(1)=min(s1,s2)/nch;
 
 sz_ch=length(res)/(prod(sz(2:end))*ni*nr);
 fprintf('Readout data points %d/%d\n',sz1_old,sz(1));
 nzero=sz_ch-sz1_old;
else
    sz1_old=sz(1)*nch;
    sz_ch=sz1_old;
    nzero=0;
    
end

nzero=nzero/2;
%res=res(1:2:end-1)+1i*res(2:2:end);
 
res=reshape(res,[2,sz_ch/2,sz(2:end)',ns,ni/ns*nr]);

res=res(1,:,:,:,:,:)+1i*res(2,:,:,:,:,:);
res=squeeze(res);
if size(res,1)~=1
 res=squeeze(res(1:sz1_old/2,:,:,:,:,:,:));
end
sz=size(res);
res=reshape(res,[sz(1)/nch,nch,sz(2:end)]);
fclose(fid);



function res=readbPar(fname,par,isnum)
%res=readbPar(fname,par,isnum)

if ~exist('isnum','var')
    isnum=true;
end

fid=fopen(fname,'r');
par=['##$',par,'='];
while 1
  b=fgetl(fid);
  if b==-1
      fclose(fid);
      error([par(4:end-1), ' not found']);
  end
  
  ind=strfind(b,par);
  if ~isempty(ind) 
    res=b(ind+length(par):end);
    break;
  end

end

if res(1)=='('
    sz=str2num(res(2:end-1));
     if numel(sz)==1
            sz=[1,sz];
     end
        
    if isnum 
     res=[];
     while 1
        
      b=fgetl(fid);
      if isnum
        tmp=str2num(b);
      else
       tmp=strread(b,'%s');  
      end
      res=[res,tmp];
      
      if length(res)==prod(sz) || ~isnum
          break;
      end
      
     end
    res=reshape(res,sz(end:-1:1));
    
    else
         b=fgetl(fid);
       res=strread(b,'%s');  
     
    end
else
    if isnum    
      res=str2double(res);
    end
    
    
end
fclose(fid);

