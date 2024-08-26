function recon_tgse_vfl_pss(fname)
%%
%{
fname='meas_MID44_tgse_vfl_EPI5_lowRes_sag_PC_FID9007.dat';
fname='meas_MID62_tgse_vfl_EPI3_lowRes_sag_PC_FID9025.dat';
fname='meas_MID58_tgse_vfl_EPI5_0_5mmInplane_TR5s_VolumeCoil_FID9120.dat';
fname = 'meas_MID18_tgse_vfl_pss_EPI3_fullFT_PE256_TR6s_FID8940.dat';
fname='meas_MID15_tgse_vfl_EPI5_lowRes_sag_lowmat_FID8959.dat';
%}
prefix=strtok(fname,'.');

lin=[];
par=[];
i=0;



while 1 
[dtmp,lintmp,partmp]=readMeasDat(fname,10000,10000*i,true);  %file name not correct

i=i+1;
%d=cat(2,d,dtmp(100:350,:));
    
if i==1
    d=dtmp;
else
    d=cat(2,d,dtmp(:,:));
end
lin=cat(2,lin,lintmp);
par=cat(2,par,partmp);


if size(dtmp,2)<10000
    break;
end

end

%%
nch=nChan_SiemensProt([prefix,'.pro']);
epi=readsPar([prefix,'.pro'],'alFree[6]');

npe_total=readsPar([prefix,'.pro'],'lPhaseEncodingLines');
pff=readsPar([prefix,'.pro'],'dSeqPhasePartialFourierForSNR');
npe=round(pff*npe_total);  % do I need to -1

pat = readsPar([prefix,'.pro'],'lAccelFact3D');

if pat==1
    lRefLines3D = 1;
else
    lRefLines3D = readsPar([prefix,'.pro'],'lRefLines3D');
end    
npe2_full=readsPar([prefix,'.pro'],'lPartitions');
npe2=lRefLines3D+npe2_full/pat-lRefLines3D/pat;

nro=size(d,1);
d=reshape(d,[nro,nch,size(d,2)/nch]);

%% check the odd-even consistency
etl=ceil(npe/epi);

npe_ref=etl*epi;
stb=readsPar([prefix,'.pro'],'lSliceTurboFactor');


%%
if pat>1
    nois_off=1;
else
    nois_off=0;
end


dref=d(:,:,1:npe_ref*stb);

lin2=lin(1:nch:end);
par2=par(1:nch:end);

lin3=lin2(npe_ref*stb+1:end);
par3=par2(npe_ref*stb+1:end);
rep=ceil(npe2/stb);



lin4=zeros(etl*epi*stb,rep);
par4=zeros(etl*epi*stb,rep);


[adc,adc2]=tgse_vfl_pss_ADCEnable(etl,epi,stb,npe,npe2);

lin4(adc>0)=lin3;
par4(adc>0)=par3;


tmp=diff(lin4,1,2);
if(any(tmp(:)>0))
    error('Phase table error');
end



dref=reshape(d(:,:,nois_off+1:nois_off+npe_ref*stb),[nro,nch,epi,stb,etl]);
data2=single(zeros(nro,nch,epi,stb,etl,rep));

data2(:,:,adc2(:)>0)=d(:,:,npe_ref*stb+1+nois_off:end);


dref(:,:,2:2:end,:,:)=flip(dref(:,:,2:2:end,:,:),1);
data2(:,:,2:2:epi-1,:,:,:)=flip(data2(:,:,2:2:epi-1,:,:,:),1);

save(['data_',prefix],'data2','dref','adc','lin4','par4');

%% plot the data and reference 
 tmp=squeeze(mean(sos(data2,2),1));

 tmp=reshape(tmp,[epi*stb*etl,rep]);
figure;imshow(tmp,[]);
figure;imshow(tmp,[0,0.0002]);
tmp=squeeze(mean(sos(dref,2),1));
figure;plot(tmp(:));
return;
%%
pc=true; %do phase correction;
use_first_pc_data=false; % use the first pc line only
corr_mag=false;


if corr_mag
    mag=abs(dref);
    mag=mean(mean(mag,1),2);
    data2=data2./repmat(mag,[nro,nch,1,1,1,rep]);
end

if pc
    an=angle(dref);
    if use_first_pc_data
         data2=data2.*exp(-1i*repmat(an(:,:,:,1),[1,1,1,stb,etl,rep]));       
    else    
         data2=data2.*exp(-1i*repmat(an,[1,1,1,1,1,rep]));   
    end
end


data3=single(zeros(nro,nch,npe,npe2));

for i=1:length(adc(:))
    if adc(i)>0
        data3(:,:,lin4(i)+1,par4(i)+1)=data2(:,:,i);
    end
end

img=ifft1c(ifft1c(data3,3),4);

img=squeeze(sos(img,2));


save_nii(make_nii(img),['recon_',prefix,'mc',num2str(corr_mag),'pc',num2str(pc),'.nii.gz']);

kdata=fft1c(data3,1);
kdata=squeeze(sos(kdata,2));

save_nii(make_nii(kdata),['kdata_',prefix,'mc',num2str(corr_mag),'pc',num2str(pc),'.nii.gz']);

%% For GRAPPA 
%{
if pat>1
    
for iro=128%1:nro

    disp(iro);
 dtmp=squeeze(ifd3(iro,:,:,:));
 
 dtmp=permute(dtmp,[2,3,1]);

 dtmp=reshape(dtmp,[npe,npe2,1,nch]);

tmp=recon_grappa2D(dtmp,par3(1,:)+1,238);
% im(:,:,iro)=tmp;
end

figure;imshow(flipud(tmp'),[0,max(tmp(:))*2/3]);
else
   
    im=ifft1c(ifft1c(ifd3,3),4);
    aim=squeeze(sos(im,2));
    save aim aim
end

%save img_recon im

%save img_150_PC Rk;
%}
%%




%%