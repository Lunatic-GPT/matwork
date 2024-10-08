function recon_FieldMap2Ds(scand)

a=readb_fid(scand);
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
a=reshape(a,[acq_size(1)/2,ns,acq_size(2),nr]);
a3=permute(a,[1,3,2,4]);

te=readbPar(fullfile(scand,'acqp'),'ACQ_echo_time');

teStep=readbPar(fullfile(scand,'method'),'teStep');

te(2)=te+teStep;

of1=readbPar(fullfile(scand,'method'),'PVM_SPackArrPhase1Offset');
of0=readbPar(fullfile(scand,'method'),'PVM_SPackArrReadOffset');
fov=readbPar(fullfile(scand,'method'),'PVM_Fov');


a=fov_shift_kspace(a3,[0,-of1],fov);


fa=fft2c(a);
fa2=flipdim(fa,1);
fa2=flipdim(fa2,2);

write_afni(abs(fa2),[scand,'_recon']);
b=angle(conj(fa2(:,:,:,1)).*fa2(:,:,:,2))/diff(te)/2/pi*1000;  %b in hertz.
write_afni(b,[scand,'_B0']);

