function recon_FLASH2D_Zspectrum(scand)

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
    chscale=[1.0,0.266,0.1879,0.1769];
else
    nch=1;
    chscale=1;
end

mtx=readbPar(fullfile(scand,'method'),'PVM_Matrix');

pe=readbPar(fullfile(scand,'acqp'),'ACQ_spatial_phase_1');
pe=pe(1:mtx(2))*mtx(2)/2;

ns=readbPar(fullfile(scand,'acqp'),'NSLICES',1);
acq_size=readbPar(fullfile(scand,'acqp'),'ACQ_size',1);
nr=readbPar(fullfile(scand,'acqp'),'NR');



a=reshape(a,[acq_size(1)/2,nch,ns,acq_size(2),nr]);

chscale=repmat(chscale,[acq_size(1)/2,1,ns,acq_size(2),nr]);
a=a.*chscale;

a2=permute(a,[1,4,3,5,2]);


ref=a2(1:end-mtx(1),:,:,:,:);
a3=a2(end-mtx(1)+1:end,:,:,:,:);
a3(:,round(pe+mtx(2)/2+1),:,:,:)=a3;


a3=phaseCorr(a3,ref(1,:,:,:,:));

of1=readbPar(fullfile(scand,'method'),'PVM_SPackArrPhase1Offset');
fov=readbPar(fullfile(scand,'method'),'PVM_Fov');
a3=fov_shift_kspace(a3,[0,-of1],fov); 


sz=size(a3);
if length(sz)<5
 sz(end+1:5)=1;
end

%a3=reshape(a3,[sz(1:3),sz(4)*sz(5)]);

fa=fft2c(a3);
fam=sqrt(mean(abs(fa(:,:,:,:,:)).^2,5));

write_afni(abs(fam),[scand,'_recon']);
%write_afni(abs(fa(:,:,:,:,2)),[scand,'_recon_ch2']);


function d2=phaseCorr(d2,nav)

  ph=angle(nav);
  
  ph=repmat(ph,[size(d2,1),1,1,1,1]);
  
 d2=d2.*exp(-1i*ph);
 
