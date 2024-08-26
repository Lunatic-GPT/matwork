function [m2,tab]=mask_cs1db_Dti(dscan)
%[m2,petable]=mask_cs1d(fid_prefix)
tab=readbPar(fullfile(dscan,'acqp'),'ACQ_spatial_phase_1');

%loadtable(tabfile(2:end-1));

Rk=readbPar(fullfile(dscan,'method'),'Rk');
sz=readbPar(fullfile(dscan,'acqp'),'ACQ_size');

nTRcs=readbPar(fullfile(dscan,'acqp'),'NR');


nDiff=readbPar(fullfile(dscan,'method'),'PVM_DwNDiffExp');

ns=readbPar(fullfile(dscan,'acqp'),'NSLICES',1);

nexpect=ns*sz(2)*nTRcs*nDiff;
if length(tab)~=nexpect
    error('length(tab)~=nexpect');
end

tab=(tab(1:nexpect)+1)*sz(2)*Rk/2+1;

tab=reshape(tab,[ns,nDiff,sz(2),nTRcs]);
tab=permute(tab,[3,1,2,4]);

tab=reshape(tab,[sz(2),ns,nTRcs*nDiff]);

mask=zeros(sz(2)*Rk,ns,nTRcs*nDiff);

for i=1:ns
    for j=1:nTRcs*nDiff
     mask(tab(:,i,j),i,j)=1;
    end
end
m2=mask;