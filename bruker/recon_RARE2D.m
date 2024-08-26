function recon_RARE2D(dname)


b=readbPar(fullfile(dname,'acqp'),'ACQ_dim');
if b~=2
    error('Not a 2D scan');
end

b=readbPar(fullfile(dname,'method'),'PVM_EncSteps1',true);
a=readb_fid(dname);
ns=readbPar(fullfile(dname,'acqp'),'NSLICES',true);
ni=readbPar(fullfile(dname,'method'),'PVM_NEchoImages');
rare=readbPar(fullfile(dname,'method'),'PVM_RareFactor');
nr=readbPar(fullfile(dname,'method'),'PVM_NRepetitions');

a2=reshape(a(:),[size(a,1),rare,ni,ns,size(a,3)/rare,nr]);
a3=permute(a2,[1,2,5,4,3,6]);

a4=reshape(a3,[size(a,1),size(a,3),ns,ni*nr]);

ra=a4;
ra(:,b-min(b)+1,:)=a4;
fra=fft2c(ra);

order=readbPar(fullfile(dname,'method'),'PVM_ObjOrderList');
fra2=fra;
fra2(:,:,order+1,:)=fra;

write_afni(abs(fra2),[dname,'_recon']);

write_afni(angle(fra2),[dname,'_reconph']);

write_afni(abs(ra),[dname,'_kspace']);



