function reorder_diffRARE2D(dname)

disp('After reorder: 1. undo recon.');
disp('               2. edit recon:');
disp('                  RECO: processing mode = FT_MODE');
disp('                  make sure RECO:output size is the same as RECO:reconstruction size...');
disp('                  RECO:quadrature options = CONJ_AND_QNEG, CONJ_AND_QNEG');
disp('               3. run reco');


if exist(fullfile(dname,'fid.orig'),'file') || exist(fullfile(dname,'acqp.orig'),'file')
  disp('File already reordered. Exit!');
  return;
end

copyfile(fullfile(dname,'fid'),fullfile(dname,'fid.orig'));

copyfile(fullfile(dname,'acqp'),fullfile(dname,'acqp.orig'));


b=readbPar(fullfile(dname,'acqp'),'ACQ_dim');
if b~=2
    error('Not a 2D scan');
end

b=readbPar(fullfile(dname,'method'),'PVM_EncSteps1',true);
[a,nzero]=readb_fid(dname);
ns=readbPar(fullfile(dname,'acqp'),'NSLICES',true);
ni=readbPar(fullfile(dname,'acqp'),'NI');
rare=readbPar(fullfile(dname,'method'),'PVM_RareFactor');

nb=ni/ns;
disp(size(a));
a2=reshape(a(:),[size(a,1),rare+2,ns,nb,size(a,2)/(rare+2)]);
a2=permute(a2,[1,2,5,3,4]);
%a3=permute(a2,[1,2,5,3,4]);

ref=a2(:,end-1:end,:,:,:);
d=a2(:,1:end-2,:,:,:);

d=phasecorr(d,ref);

fmt=readbPar(fullfile(dname,'acqp'),'GO_raw_data_format',false);

if ~strcmp(fmt,'GO_32BIT_SGN_INT')
  error('fid data format %s not supported',fmt);
end

d=reshape(d,[size(d,1),size(d,2)*size(d,3),ns,nb]);
ra=d;
ra(:,b-min(b)+1,:,:)=d;

%ra=flipdim(ra,1);
ra(end+1:end+nzero,:,:)=0;
%ra(end+1:end+51,end+1:end+24,:,:)=0;
fra=fft2c(ra);
figure;imshow(abs(fra(:,:,1,1,1)),[]);
ra=permute(ra,[1,3,2,4]);

fid=fopen(fullfile(dname,'fid'),'w','ieee-le');
ra=flipdim(ra,3);
d=ra(:);
d2=cat(2,real(d),imag(d));
d2=d2';

disp(size(ra));

fwrite(fid,d2(:),'int32');

changebPar(fullfile(dname,'acqp'),'NR',1,nb);
changebPar(fullfile(dname,'acqp'),'NI',1,ns);
changebPar(fullfile(dname,'acqp'),'ACQ_size',1,[(size(ra,1)-nzero)*2,size(ra,3)]);
changebPar(fullfile(dname,'acqp'),'ACQ_phase_factor',1,1);
changebPar(fullfile(dname,'acqp'),'ACQ_phase_factor',1,1);

changebPar(fullfile(dname,'acqp'),'ACQ_spatial_size_1',1,size(ra,2));
od=readbPar(fullfile(dname,'acqp'),'ACQ_obj_order');
changebPar(fullfile(dname,'acqp'),'ACQ_obj_order',1,od(1:8));
changebPar(fullfile(dname,'acqp'),'ACQ_phase_encoding_mode',0,{'Read','Linear'});

fclose(fid);


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














