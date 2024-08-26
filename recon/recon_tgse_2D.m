function recon_tgse_2D(fname)
%%

%fname = 'meas_MID68_tgse_FID9130.dat';
%fname = 'meas_MID69_tgse_NA1_FID9131.dat';

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


if size(dtmp,2)<10000
    break;
end

end
%%
ns=readsPar([prefix,'.pro'],'sSliceArray.lSize');

epi=readsPar([prefix,'.pro'],'lEPIFactor');
tb=readsPar([prefix,'.pro'],'lTurboFactor');
na = readsPar([prefix,'.pro'],'lAverages');  

lin=reshape(lin,[epi,tb,ns,length(lin(:))/ns/epi/tb]);

nro=size(d,1);
nrep=size(lin,4);
d=reshape(d,[nro,epi,tb,ns,nrep]);

fd=fft1c(d,1);
fd(:,2:2:epi-1,:,:,:,:)=flip(fd(:,2:2:epi-1,:,:,:,:),1);

ifd=ifft1c(fd,1);

d2=ifd(:,:,:,:,2:end);
dref=ifd(:,:,:,:,1);
lin2=reshape(lin(:,:,1,2:end),[epi*tb,na,(nrep-1)/na]);

d2=reshape(d2,[nro,epi*tb,ns,na,(nrep-1)/na]);

d2=permute(d2,[1,2,4,5,3]);

im=zeros(nro,epi*tb*(nrep-1)/na,ns,na);
for i=1:na
    
    tmp=reshape(d2(:,:,i,:,:),[nro,epi*tb*(nrep-1)/na,ns]);
    lintmp=lin2(:,i,:);
    
    tmp(:,lintmp+1,:)=tmp;
    
    im(:,:,:,i)=ifft1c(tmp,2);
end

%%
dref_max=max(abs(dref),[],1);

mag=mean(abs(dref),1);
figure;plot(abs(dref(:,:,1,6))./repmat(dref_max(:,:,1,6),[512,1,1,1]));








