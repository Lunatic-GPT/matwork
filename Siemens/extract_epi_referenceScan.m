function extract_epi_referenceScan(fname)

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

if i==8
    break;
end

if size(dtmp,2)<10000
    break;
end

end
%%
nslice=readsPar([prefix,'.pro'],'sSliceArray.lSize');

nch=nChan_SiemensProt([prefix,'.pro']);


lAccelFactPE=readsPar([prefix,'.pro'],'lAccelFactPE');

if lAccelFactPE >1
    paton=1;
else
    paton=0;
end

npe_ref=3;
npe=31;



if paton
    fd=fft1c(d(:,nch+1:npe*nch*nslice+nch),1);
    
else
    fd=fft1c(d(:,1:npe*nch*nslice),1);
    
end
    
fd=reshape(fd,[size(d,1),nch,npe,nslice]);

fd(:,:,2:2:end,:)=flip(fd(:,:,2:2:end,:),1);
ifd=ifft1c(fd,1);

ifd=ifd(:,:,1:npe_ref,:);

figure;imshow(abs(squeeze(ifd(:,1,:,:))),[]);

stotal=mean(sos(ifd,2),1);
figure;plot(squeeze(stotal));
prefix = strtok(fname,'.');

save(['ref_data_',prefix], 'ifd');
