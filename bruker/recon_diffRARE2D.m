function recon_diffRARE2D(dname)


b=readbPar(fullfile(dname,'acqp'),'ACQ_dim');
if b~=2
    error('Not a 2D scan');
end

b=readbPar(fullfile(dname,'method'),'PVM_EncSteps1',true);
a=readb_fid(dname);
ns=readbPar(fullfile(dname,'acqp'),'NSLICES',true);
ni=readbPar(fullfile(dname,'acqp'),'NI');
rare=readbPar(fullfile(dname,'method'),'PVM_RareFactor');

nb=ni/ns;

a2=reshape(a(:),[size(a,1),rare+2,ns,nb,size(a,2)/(rare+2)]);
a3=permute(a2,[1,2,5,3,4]);

ref=a3(:,end-1:end,:,:,:);
d=a3(:,1:end-2,:,:,:);

d=phasecorr(d,ref);
d=reshape(d,[size(d,1),size(d,2)*size(d,3),ns,nb]);

ra=d;
ra(:,b-min(b)+1,:,:)=d;
ra=flipdim(ra,1);
fra=fft2c(ra);

order=readbPar(fullfile(dname,'method'),'PVM_ObjOrderList');
fra2=fra;
fra2(:,:,order+1,:)=fra;

write_afni(abs(fra2),[dname,'_recon2']);

function res=phasecorr(d,ref)

if mod(size(d,2),2)~=0
    error('Only works for even rare factors');
end
res=zeros(size(d));

fd=fft1c(d,1);
fref=fft1c(ref,1);

for i=1:2
  ph=angle(fref(:,i,:,:,:));
  ph=repmat(ph,[1,size(res,2)/2,1,1,1]);
  
  res(:,i:2:end,:,:,:)=fd(:,i:2:end,:,:,:).*exp(-1i*ph);
end

res=ifft1c(res,1);














