function extract_tgse_referenceScan(fname)

prefix=strtok(fname,'.');

d=[];
lin=[];
par=[];
i=0;

while 1 
[dtmp,lintmp,partmp]=readMeasDat(fname,10000,10000*i,true);  %file name not correct

i=i+1;
%d=cat(2,d,dtmp(100:350,:));
d=cat(2,d,dtmp(:,:));
lin=cat(2,lin,lintmp);
par=cat(2,par,partmp);

if i==2
    break;
end

if size(dtmp,2)<10000
    break;
end

end
%%

nch=nChan_SiemensProt([prefix,'.pro']);

epi=readsPar([prefix,'.pro'],'alFree[6]');
if isempty(epi)
    epi=1;
end

lAccelFactPE=readsPar([prefix,'.pro'],'lAccelFactPE');

if lAccelFactPE >1
    paton=1;
else
    paton=0;
end

npe_total=readsPar([prefix,'.pro'],'lPhaseEncodingLines');
pff=readsPar([prefix,'.pro'],'dSeqPhasePartialFourierForSNR');
te=readsPar([prefix,'.pro'],'alTE[0]');
npe=round(pff*npe_total);
%npe=182;
npe_ref=ceil(npe/epi)*epi;

echoSpacing=te/(npe-(npe_total-1)/2);  %only works for epi=1

if paton
    fd=fft1c(d(:,nch+1:npe_ref*nch+nch),1);
    
    fd_highk=fft1c(d(:,npe_ref*nch+nch+1:end),1);

else
    fd=fft1c(d(:,1:npe_ref*nch),1);
    
    fd_highk=fft1c(d(:,npe_ref*nch+1:end),1);
end
    
fd=reshape(fd,[size(d,1),nch,epi,npe_ref/epi]);

fd(:,:,2:2:epi-1,:)=flip(fd(:,:,2:2:epi-1,:),1);
ifd=ifft1c(fd,1);

ifd=reshape(ifd,[size(d,1),nch,npe_ref]);  %

%figure;imshow(abs(squeeze(ifd(:,1,:))),[]);

ich=1;

figure;plot(vec(abs(ifd(:,ich,:))));
prefix = strtok(fname,'.');

save(['ref_xdata_',prefix], 'ifd');


