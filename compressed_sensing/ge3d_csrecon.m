function ge3d_csrecon(fid_prefix,reg)

if ~exist('reg','var')
 reg=0.05;
end
% this is not yet tested;
tabfile=readPar(fid_prefix,'petable');

if exist(tabfile(2:end-1),'file')
    tab=load(tabfile(2:end-1));
elseif strcmp(computer,'PCWIN64')
   tab=load(['z:/tmp.home/xiaopeng/newvnmrsys/tablib/',tabfile(2:end-1)]);
else
   tab=load(['/home/xiaopeng/vnmrsys/tablib/',tabfile(2:end-1)]);
end

cs_pe=readPar(fid_prefix,'cs_pe');
cs_pe2=readPar(fid_prefix,'cs_pe2');
np=readPar(fid_prefix,'np');
nv=readPar(fid_prefix,'nv');
nv2=readPar(fid_prefix,'nv2');
arraydim= readPar(fid_prefix,'arraydim');
nTRcs=readPar(fid_prefix,'nTRcs');
tab2=zeros(cs_pe*cs_pe2,arraydim);

for i=1:arraydim  
  os=mod((i-1),nTRcs);  
  tab2(:,i)=tab(os*cs_pe*cs_pe2+1:(os+1)*cs_pe*cs_pe2);  
end


mask=zeros(nv,nv2,arraydim);

for i=1:arraydim
    tmp=zeros(nv,nv2);
    tmp(tab2(:,i))=1;
  mask(:,:,i)=tmp;
end
m2=shiftdim(mask,-1);
m2=repmat(m2,[np/2,1,1,1]);

z=read_fid([fid_prefix,'.fid']);
z2=zeros(np/2,nv,nv2,arraydim);
z2(m2>0)=z;
z2=z2/max(abs(z2(:)));
XFM=FTt;
z_init=mean_sampled_kt(z2,m2>0);
z_init=repmat(z_init,[1,1,1,size(z2,4)]);
z_init=z2;

fz2=ifft1c(z2,1);
im_res=zeros(size(fz2));

for i=1:size(fz2,1)
 ztmp=reshape(fz2(i,:,:,:),[size(fz2,2),size(fz2,3),1,size(fz2,4)]);
 mtmp=reshape(m2(i,:,:,:),[size(fz2,2),size(fz2,3),1,size(fz2,4)]);
 im_res(i,:,:,:)=run_cs(ztmp,ztmp,[],[],mtmp,XFM,0,reg,15,true);
end

write_afni(abs(im_res),sprintf('%s_reg%4.3f',fid_prefix,reg));

write_afni(angle(im_res),sprintf('%s_reg%4.3f_ph',fid_prefix,reg));




