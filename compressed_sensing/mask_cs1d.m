function [m2,tab]=mask_cs1d(fid_prefix)
%[m2,petable]=mask_cs1d(fid_prefix)
tabfile=readPar(fid_prefix,'petable');
if ~exist(fullfile(pwd,tabfile(2:end-1)),'file') 
  fprintf('use %s saved in vnmrsys/tablib\n',tabfile(2:end-1)); 
  tab=load(['~/vnmrsys/tablib/',tabfile(2:end-1)]);
else
  tab=load(fullfile(pwd,tabfile(2:end-1)));  
end

%loadtable(tabfile(2:end-1));

cs_pe=readPar(fid_prefix,'npecs');
np=readPar(fid_prefix,'np');
nv=readPar(fid_prefix,'nv');
nTRcs=readPar(fid_prefix,'nTRcs');
ne=readPar(fid_prefix,'arraydim');
ns=readPar(fid_prefix,'ns');

if length(tab)~=ns*cs_pe*nTRcs
    warning('length(tab)~=ns*cs_pe*nTRcs');
end

tab=tab(1:cs_pe*ns*nTRcs)+nv/2+1;

tab=reshape(tab,[cs_pe,ns,nTRcs]);

mask=zeros(nv,ns,nTRcs);

for i=1:ns
    for j=1:ne
        iTRcs=mod(j-1,nTRcs)+1;
     mask(tab(:,i,iTRcs),i,j)=1;
    end
end
m2=mask;
%mask=shiftdim(mask,-1);
%m2=repmat(mask,[np/2,1,1,1]);