function [res,nzero]=readb_fid(d)


if ~exist(d,'dir')  % d is a file
    [d,fname,suf]=fileparts(d);
    fname=[fname,suf];
else
    fname='fid';
end

sz=readbPar(fullfile(d,'acqp'),'ACQ_size',true);
nr=readbPar(fullfile(d,'acqp'),'NR',true);
ni=readbPar(fullfile(d,'acqp'),'NI',true);

if numel(sz)>1
 ns = readbPar(fullfile(d,'acqp'),'NSLICES',true);
else
    ns=1;
end

fid=fopen(fullfile(d,fname),'r','ieee-le');

fmt=readbPar(fullfile(d,'acqp'),'GO_raw_data_format',false);

if strcmp(fmt,'GO_32BIT_SGN_INT')
 res=fread(fid,'int32');
 
else
    fprintf('%s :',fmt);
    fclose(fid);
    error('unknown format');
end

acqmod=readbPar(fullfile(d,'acqp'),'ACQ_experiment_mode',false);
if ~strcmp('SingleExperiment',acqmod)
    rcvrs=readbPar(fullfile(d,'acqp'),'ACQ_ReceiverSelect',false);
    nyes=strmatch('Yes',rcvrs);
    nch=length(nyes);    
   
else
    nch=1;
end

bs=readbPar(fullfile(d,'acqp'),'GO_block_size',false);

if strcmp('Standard_KBlock_Format',bs);
    sz1_old=sz(1)*nch;
 
 %s1=2.^ceil(log2(sz(1)*nch));
 %s2=2.^ceil(log2(sz(1)*nch/10))*10;
 %sz(1)=min(s1,s2)/nch;
 
 sz_ch=length(res)/(prod(sz(2:end))*ni*nr);
 fprintf('Readout data points %d/%d\n',sz1_old,sz(1));
 nzero=sz_ch-sz1_old;
else
    sz1_old=sz(1)*nch;
    sz_ch=sz1_old;
    nzero=0;
    
end

nzero=nzero/2;
%res=res(1:2:end-1)+1i*res(2:2:end);
 
res=reshape(res,[2,sz_ch/2,sz(2:end)',ns,ni/ns*nr]);

res=res(1,:,:,:,:,:)+1i*res(2,:,:,:,:,:);
res=squeeze(res);
if size(res,1)~=1
 res=squeeze(res(1:sz1_old/2,:,:,:,:,:,:));
end
sz=size(res);
res=reshape(res,[sz(1)/nch,nch,sz(2:end)]);
fclose(fid);




