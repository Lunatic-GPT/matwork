function recon_DtiStandard(scand)
%gems_csrecon(fid_prefix,method[,ref,ref_ind,recon_ind])

z=readb_fid(scand);


b=readbPar(fullfile(scand,'acqp'),'ACQ_dim');
if b~=2
    error('Not a 2D scan');
end

acqmod=readbPar(fullfile(scand,'acqp'),'ACQ_experiment_mode',false);
if ~strcmp('SingleExperiment',acqmod)
    rcvrs=readbPar(fullfile(scand,'acqp'),'ACQ_ReceiverSelect',false);
    nyes=strmatch('Yes',rcvrs);
    nch=length(nyes);    
   
else
    nch=1;
end

ns=readbPar(fullfile(scand,'acqp'),'NSLICES',1);
acq_size=readbPar(fullfile(scand,'acqp'),'ACQ_size',1);
nr=readbPar(fullfile(scand,'acqp'),'NR');
nq=readbPar(fullfile(scand,'method'),'PVM_DwNDiffExp');

of1=readbPar(fullfile(scand,'method'),'PVM_SPackArrPhase1Offset');

of0=readbPar(fullfile(scand,'method'),'PVM_SPackArrReadOffset');

fov=readbPar(fullfile(scand,'method'),'PVM_Fov');

z=reshape(z,[acq_size(1)/2,nch,ns,nq,acq_size(2),nr]);

z=permute(z,[1,5,3,2,4,6]);

z=mean(z,6);

sz=size(z);

z=reshape(z,[sz(1:3),prod(sz(4:end))]);

z=fov_shift_kspace(z,[-of0,-of1],fov);
fz=ifft2c(z);

write_afni(abs(fz),[scand,'_mag']);

write_afni(angle(fz),[scand,'_ph']);



