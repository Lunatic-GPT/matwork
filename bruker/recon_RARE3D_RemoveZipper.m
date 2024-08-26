function recon_RARE3D_RemoveZipper(dname)
% recon_RARE3D_RemoveZipper dname
% dname: the scan data folder (usually a number)
% Example: recon_RARE3D_RemoveZipper 23
% The output image is in .mat format

b=readbPar(fullfile(dname,'acqp'),'ACQ_dim');
if b~=3
    error('Not a 3D scan');
end

b=readbPar(fullfile(dname,'method'),'PVM_EncSteps1',true);
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
%a=reshape(a,[acq_size(1)/2,nch,acq_size(2),acq_size(3),nr*ni]);
%a2=permute(a,[1,3,4,2,5]);

a=reshape(a,[acq_size(1)/2,nch,acq_size(2),ni*nr,acq_size(3)]);
a2=permute(a,[1,3,5,4,2]);
sz=size(a2);
if length(sz)<5
 sz(end+1:5)=1;
end
a2=reshape(a2,[sz(1:3),sz(4)*sz(5)]);

ra=a2;
ra(:,b-min(b)+1,:,:)=a2;


fra2=fft1c(ra,2);
fra23=fft1c(fra2,3);

sz=size(ra);
fra23([1,end],sz(2)/2+1,sz(3)/2+1)=0;
fra123=fft1c(fra23,1);

afra123=abs(fra123);

save([dname,'_recon'],'afra123');








