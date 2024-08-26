function recon_FieldMap2D(scand)

for i=1:2
a=readb_fid(scand{i});
b=readbPar(fullfile(scand{i},'acqp'),'ACQ_dim');
if b~=2
    error('Not a 2D scan');
end

acqmod=readbPar(fullfile(scand{i},'acqp'),'ACQ_experiment_mode',false);
if ~strcmp('SingleExperiment',acqmod)
    rcvrs=readbPar(fullfile(scand{i},'acqp'),'ACQ_ReceiverSelect',false);
    nyes=strmatch('Yes',rcvrs);
    nch=length(nyes);    
   
else
    nch=1;
end

ns=readbPar(fullfile(scand{i},'acqp'),'NSLICES',1);
acq_size=readbPar(fullfile(scand{i},'acqp'),'ACQ_size',1);
nr=readbPar(fullfile(scand{i},'acqp'),'NR');
a=reshape(a,[acq_size(1)/2,nch,ns,acq_size(2),nr]);

a2=permute(a,[1,4,3,2,5]);
sz=size(a2);
if length(sz)<5
 sz(end+1:5)=1;
end
a2=reshape(a2,[sz(1:3),sz(4)*sz(5)]);
a2=mean(a2,4);
a3(:,:,:,i)=a2;
end  


for i=1:2
te(i)=readbPar(fullfile(scand{i},'acqp'),'ACQ_echo_time');

end


of1=readbPar(fullfile(scand{1},'method'),'PVM_SPackArrPhase1Offset');
of0=readbPar(fullfile(scand{1},'method'),'PVM_SPackArrReadOffset');
fov=readbPar(fullfile(scand{1},'method'),'PVM_Fov');


a=fov_shift_kspace(a3,[0,-of1],fov);


fa=fft2c(a);
fa2=flipdim(fa,1);
fa2=flipdim(fa2,2);

write_afni(abs(fa2),[scand{1},'_recon']);
b=angle(conj(fa2(:,:,:,1)).*fa2(:,:,:,2))/diff(te)/2/pi*1000;  %b in hertz.
write_afni(b,[scand{1},'_B0']);

