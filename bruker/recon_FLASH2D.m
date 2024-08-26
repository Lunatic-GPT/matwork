function recon_FLASH2D(scand)

a=readb_fid(scand);
b=readbPar(fullfile(scand,'acqp'),'ACQ_dim');
if b~=2
    error('Not a 2D scan');
end

%sz=size(a);
%a=reshape(a,[sz(1),sz(3),sz(2),sz(4)]);
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
a=reshape(a,[acq_size(1)/2,nch,ns,acq_size(2),nr]);

a2=permute(a,[1,4,3,2,5]);
sz=size(a2);
if length(sz)<5
 sz(end+1:5)=1;
end
a2=reshape(a2,[sz(1:3),sz(4)*sz(5)]);

fa=fft2c(a2);
%fa2=fft1c(fa,3);

write_afni(abs(fa),[scand,'_recon']);
