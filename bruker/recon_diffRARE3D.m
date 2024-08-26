function recon_diffRARE3D(dname)


b=readbPar(fullfile(dname,'acqp'),'ACQ_dim');
if b~=3
    error('Not a 3D scan');
end

b=readbPar(fullfile(dname,'method'),'PVM_EncSteps1',true);

b2=readbPar(fullfile(dname,'method'),'PVM_EncSteps2',true);

a=readb_fid(dname);
%rare=readbPar(fullfile(dname,'method'),'PVM_RareFactor');


%a2=reshape(a(:),[size(a,1),rare,size(a,2)/rare,size(a,3),size(a,4)]);

%a4=reshape(a3,[size(a,1),size(a,2),ns,ni]);

acqmod=readbPar(fullfile(dname,'acqp'),'ACQ_experiment_mode',false);
if ~strcmp('SingleExperiment',acqmod)
    rcvrs=readbPar(fullfile(dname,'acqp'),'ACQ_ReceiverSelect',false);
    nyes=strmatch('Yes',rcvrs);
    nch=length(nyes);    
   
else
    nch=1;
end

acq_size=readbPar(fullfile(dname,'acqp'),'ACQ_size',1);
nr=readbPar(fullfile(dname,'acqp'),'NR');
ni=readbPar(fullfile(dname,'acqp'),'NI');
rare=readbPar(fullfile(dname,'acqp'),'ACQ_rare_factor');

%{
a=reshape(a,[acq_size(1)/2,nch,rare,ni*nr,acq_size(2)/rare,acq_size(3)]);
a2=permute(a,[1,3,5,6,4,2]);
a2=reshape(a2,[acq_size(1)/2,acq_size(2),acq_size(3),ni*nr*nch]);
%}
a2=reshape(a,[acq_size(1)/2,acq_size(2),acq_size(3),ni*nr]);


ra=a2;
ra(:,b-min(b)+1,b2-min(b2)+1,:)=a2;


fra=fft2c(ra);

fra2=fft1c(fra,3);
%fra2=permute(fra2,[2,3,1,4]);

%ra=permute(ra,[2,3,1,4]);

write_afni(abs(fra2),[dname,'_recon']);
write_afni(abs(ra),[dname,'_kspace']);




