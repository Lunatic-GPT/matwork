function array_epics(pss,nTR,outf,off_on_off)
%array_epics(pss,nTR,outf)
% seqcon = 'nsnnn'
% nTR excludes the first reference scan.

ns=length(pss);
image=ones(ns,nTR+1);

image(:,1)=0;


pss=repmat(pss(:),[1,nTR+1]);


fid=fopen(outf,'w');

for i=1:length(pss(:))
 fprintf(fid,'%f %d\n',pss(i),image(i));
end

fclose(fid);

if ~exist('off_on_off','var')
  return;
end

sli_val=zeros(ns,nTR+1);
ntrl=nTR/sum(off_on_off);

if mod(nTR,sum(off_on_off))~=0
    error('mod(nTR,sum(off_on_off))~=0');
end

for i=1:ntrl
     sli_val(:,(i-1)*sum(off_on_off)+off_on_off(1)+2:(i-1)*sum(off_on_off)+sum(off_on_off(1:2))+1)=1;      
end

fid=fopen([outf,'_sli_val'],'w');
for i=1:length(pss(:))
 fprintf(fid,'%d %d\n',sli_val(i),image(i));
end

fclose(fid);


