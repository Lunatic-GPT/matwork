function res = recon_EPImb_Siemens(filename,varargin)
% Example:
%
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70]);
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70],'Echo',2,'Rep',50);
% this is tested for the MGH sequence; assuming single-band data comes
% before the grappa calibration data.

tic;
prefix=strtok(filename,'.');

[a,kindex]=readMeasDat(filename,'e:/dropbox/matwork/Siemens/coilScaling.txt',[]);

kindex=kindex(1:32:end)+1;
a=ifft1c(a,1);
a=a(end/4+1:end-end/4,:);
a=fft1c(a,1);
%%
    % MDH Header Information Retrieval Ends Here...
   
    nro=size(a,1);
    npe=readsPar([prefix,'.pro'],'lPhaseEncodingLines');
    
    nb=readsPar([prefix,'.pro'],'adFree[1]');  % maybe adFree[2]
    % Slices    

 
nch=readsPar([prefix,'.pro'],'lMaximumNofRxReceiverChannels');

    nsl=readsPar([prefix,'.pro'],'asSlice');
    nsl=length(nsl)/7;
    
    nrep=readsPar([prefix,'.pro'],'lRepetitions');
    nref=readsPar([prefix,'.pro'],'lRefLinesPE');
    accel=readsPar([prefix,'.pro'],'lAccelFactPE');
    nrep=nrep+1;
   
     a=reshape(a,[nro,nch,size(a,2)/nch]);
 a=permute(a,[1,3,2]);
 
 %% noise
 dnois=a(:,1,:);
 
 
   %% single band
    k_center=kindex(2);
    k_center_ind=find(kindex==k_center);
    
    l1b=(k_center_ind(5)-k_center_ind(3)+2)*nsl;
    d1b=reshape(a(:,2:1+l1b,:),[nro,l1b/nsl,nsl,nch]);
    
    if mod(nsl,2)==0
     sorder=[2:2:nsl,1:2:nsl-1];
    else
     sorder=[1:2:nsl-1,2:2:nsl];
    end
    d1b(:,:,sorder,:)=d1b;
  
    
    d1b=epi_flipro_phasecor(d1b);
    k1b=kindex(5:4+l1b/nsl-3);
    
    %% GRAPPA calibration
    lref=(nref+3*accel)*nsl;
    dref=reshape(a(:,2+l1b:lref+1+l1b,:),[nro,nref/accel+3,nsl,accel,nch]);
    dref=permute(dref,[1,2,4,3,5]);
  
    
    dref=epi_flipro_phasecor(dref);
   
    
    dref=reshape(dref,[nro,nref,nsl,nch]);
    dref(:,:,sorder,:,:)=dref;
   
    kref=kindex(2+l1b:lref+1+l1b);
    kref=reshape(kref,[nref/accel+3,nsl,accel]);
    kref=kref(4:end,1,:);
    kref=kref(:);
    kref=kref-min(kref)+1;

    dref(:,kref,:,:)=dref;
   %%
   
d1b_2=zeros(nro,max(k1b)+accel-1,nsl,nch);
d1b_2(:,k1b,:,:)=d1b;

Rkall=recon_grappa_1rep(d1b_2,dref,k1b);


Rkall=Rkall(:,1:end-2,:,:);
[img,ph]=partialFT(Rkall,k_center-1);
img1b=sqrt(sum(abs(img).^2,4));

save(['recon_',prefix,'_1b.mat'],'img1b');
%% multi-band     
   lmb=l1b/nb*nrep;
   
    dmb=reshape(a(:,2+l1b+lref:1+l1b+lref+lmb,:),[nro,lmb*nb/nsl/nrep,nsl/nb,nrep,nch]);
    
    dmb=epi_flipro_phasecor(dmb);
     
     
    if mod(nsl/nb,2)==0
     sorder=[2:2:nsl/nb,1:2:nsl/nb-1];
    else
     sorder=[1:2:nsl/nb+1,2:2:nsl/nb];
    end
    dmb(:,:,sorder,:,:)=dmb;
   
   
dmb_2=zeros(nro,max(k1b)+accel-1,nsl/nb,nrep,nch);
dmb_2(:,k1b,:,:,:)=dmb;
    
    for i=1:size(dmb_2,3)
    
        x=reshape(dmb_2(:,:,2,1,:),[size(dmb_2,1),size(dmb_2,2),1,size(dmb_2,5)]);
        Rkall=recon_grappa_1rep(x,sum(dref(:,:,[2,5],:),3),k1b);
        
        
    end
    
    
method='slice_grappa';
 switch method
     
     case 'slice_grappa'
         
    kernY=-2:2;
    kernX=-2:2;
    
     dmb4=zeros(nro,npe,nsl,nch);
    for i=1:nsl/nb 
        for j=1:nb
         coef=myGRAPPA_2data_xp(squeeze(sum(d1b(:,:,i:end/nb:end,:),3)),squeeze(d1b(:,:,i+(j-1)*nsl/nb,:)),1:size(d1b,2),1:size(d1b,2),kernY,kernX,0,true);  
         dmb_sep=myGRAPPA_2data_xp(squeeze(dmb(:,:,i,:)),[],1:size(dmb,2),1:size(dmb,2),kernY,kernX,coef);
         dmb_sep=reshape(dmb_sep,[size(dmb_sep,1),size(dmb_sep,2),1,size(dmb_sep,3)]); 
         dmb_sep2=phase_corr(dmb_sep,ref1b(:,:,i+(j-1)*nsl/nb,:));
         isamp=2:3:116;  
         dmb4(:,isamp,i+(j-1)*nsl/nb,:)=flipdim(dmb_sep2,1);   
         recon_grappa(dmb4(:,:,i+(j-1)*nsl/nb,:),dref3(:,:,i+(j-1)*nsl/nb,:),isamp);
    end
    end
    
     case 'cat_slice'
    for ir=1
        
    dmb2=squeeze(dmb(:,:,:,:,ir));
                                                   % multi-band data should
                                                 % be nro*39*nch
    fd1b=fft(d1b,[],3);
    fd1b=circshift(fd1b,[0,0,ceil(nsl/2),0]);
    
    dmb3=repmat(zeros(size(dmb2)),[1,1,4,1]);
    ismp=1:4:size(dmb3,3);
    dmb3(:,:,1:4:end,:)=fft(dmb2,[],3);
    
    kernY=[-5,-1,3,7];
    kernX=-2:2;
   
    len=size(d1b,2);
    
    % tic; [tmp,tmp2,coef]=myGRAPPA_xp(squeeze(fd1b(:,round(1),:,:)),1:nsl,1:nsl,kernY,kernX,0,true);toc;
    for k=1:len
     [tmp,tmp2,coef]=myGRAPPA_xp(squeeze(fd1b(:,k,:,:)),1:nsl,1:nsl,kernY,kernX,0);
      Rk(:,k,:,:)=myGRAPPA_xp(squeeze(dmb3(:,k,:,:)),ismp,1:nsl,kernY,kernX,coef);
    end
    
    kernY=[-6,-2,2,6];
    
    %  [tmp,tmp2,coef]=myGRAPPA_xp(squeeze(fd1b(:,round(end/2),:,:)),1:nsl,1:nsl,kernY,kernX,0);
    for k=1:len
      [tmp,tmp2,coef]=myGRAPPA_xp(squeeze(fd1b(:,k,:,:)),1:nsl,1:nsl,kernY,kernX,0);
      Rk2(:,k,:,:)=myGRAPPA_xp(squeeze(dmb3(:,k,:,:)),ismp,1:nsl,kernY,kernX,coef);
    end
    Rk(:,:,ismp+2,:)=Rk2(:,:,ismp+2,:);
    
    
    kernY=[-7,-3,1,5];
    
     % [tmp,tmp2,coef]=myGRAPPA_xp(squeeze(fd1b(:,round(end/2),:,:)),1:nsl,1:nsl,kernY,kernX,0);
    for k=1:len
      [tmp,tmp2,coef]=myGRAPPA_xp(squeeze(fd1b(:,k,:,:)),1:nsl,1:nsl,kernY,kernX,0);
      Rk3(:,k,:,:)=myGRAPPA_xp(squeeze(dmb3(:,k,:,:)),ismp,1:nsl,kernY,kernX,coef);
    end
    Rk(:,:,ismp+3,:)=Rk3(:,:,ismp+3,:);
    
    Rk=ifft(Rk,[],3);
    
    %dmb_sep=permute(dmb_sep,[1,3,2,4]);
%dmb_sep=Rk;        
      dmb_sep=phase_corr(Rk,ref1b);
      
     dmb4=zeros(nro,npe,nsl,nch);

     isamp=2:3:116;
        
      dmb4(:,isamp,:,:)=dmb_sep;   
     
      dmb4=flipdim(dmb4,1);
      recon_grappa(dmb4(:,:,8,:),dref3(:,:,8,:),isamp);
 
    end
     otherwise
        
 end
         
 
 
  
    
 
fprintf('Read data took %s s\n',toc);

function d1bk=phase_corr(d1b,ref)
    d1b=ifft1c(d1b,1);
    ref=ifft1c(ref,1);
    ph=angle(ref);
 
    d1b(:,1:2:end,:,:,:)=d1b(:,1:2:end,:,:,:).*repmat(exp(-1i*ph(:,1,:,:,:)),[1,ceil(size(d1b,2)/2),1,1,1]);
    d1b(:,2:2:end,:,:,:)=d1b(:,2:2:end,:,:,:).*repmat(exp(-1i*ph(:,2,:,:,:)),[1,floor(size(d1b,2)/2),1,1,1]);
     d1bk=fft1c(d1b,1);
%}
function Rkall=recon_grappa_1rep(data,dref,isamp)


 nro=size(data,1);
 npe=size(data,2);
 nsl=size(data,3);
 nch=size(data,4);
 
 Rkall=zeros(nro,npe,nsl,nch);
  
  kernelX=[-2,-1,0,1,2];
 kernY=[-4,-1,2,5];
   
    
     
    for irep=1:nsl
    disp(irep);
     [tmp,coef1]=myGRAPPA_xp(squeeze(dref(:,:,irep,:)),1:size(dref,2),1:size(dref,2),kernY,kernelX,0);
    
    
     kernY=[-5,-2,1,4];
     [tmp,coef2]=myGRAPPA_xp(squeeze(dref(:,:,irep,:)),1:size(dref,2),1:size(dref,2),kernY,kernelX,0);
     
    %    data=fft1c(squeeze(dref3(:,:,isl,:)),1);
   
    kernY=[-4,-1,2,5];
    
    Rk=myGRAPPA_xp(squeeze(data(:,:,irep,:)),isamp,0,kernY,kernelX,coef1);
    
     kernY=[-5,-2,1,4];
    
     Rk2=myGRAPPA_xp(squeeze(data(:,:,irep,:)),isamp,0,kernY,kernelX,coef2);
    
     Rk(:,isamp+2,:)=Rk2(:,isamp+2,:);
     Rkall(:,:,irep,:)=Rk;
    end
%{
function mdhBitFields = determineBitFields(evalInfo)

bits = num2str(dec2bin(evalInfo));
bits = fliplr([repmat('0',1,32-length(bits)),bits]);
setFlags = find(bits=='1')-1;

mdhBitFields.MDH_ACQEND = any(setFlags==0);
mdhBitFields.MDH_RTFEEDBACK = any(setFlags==1);
mdhBitFields.MDH_HPFEEDBACK = any(setFlags==2);
mdhBitFields.MDH_ONLINE    = any(setFlags==3);
mdhBitFields.MDH_OFFLINE   = any(setFlags==4);
mdhBitFields.MDH_LASTSCANINCONCAT = any(setFlags==8);       % Flag for last scan in concatination
mdhBitFields.MDH_RAWDATACORRECTION = any(setFlags==10);      % Correct the rawadata with the rawdata correction factor
mdhBitFields.MDH_LASTSCANINMEAS = any(setFlags==11);      % Flag for last scan in measurement
mdhBitFields.MDH_SCANSCALEFACTOR = any(setFlags==12);      % Flag for scan specific additional scale factor
mdhBitFields.MDH_2NDHADAMARPULSE = any(setFlags==13);      % 2nd RF exitation of HADAMAR
mdhBitFields.MDH_REFPHASESTABSCAN = any(setFlags==14);      % reference phase stabilization scan
mdhBitFields.MDH_PHASESTABSCAN = any(setFlags==15);      % phase stabilization scan
mdhBitFields.MDH_D3FFT     = any(setFlags==16);      % execute 3D FFT
mdhBitFields.MDH_SIGNREV   = any(setFlags==17);      % sign reversal
mdhBitFields.MDH_PHASEFFT  = any(setFlags==18);      % execute phase fft
mdhBitFields.MDH_SWAPPED   = any(setFlags==19);      % swapped phase/readout direction
mdhBitFields.MDH_POSTSHAREDLINE = any(setFlags==20);      % shared line
mdhBitFields.MDH_PHASCOR   = any(setFlags==21);      % phase correction data
mdhBitFields.MDH_PATREFSCAN = any(setFlags==22);      % additonal scan for PAT reference line/partition
mdhBitFields.MDH_PATREFANDIMASCAN = any(setFlags==23);      % additonal scan for PAT reference line/partition that is also used as image scan
mdhBitFields.MDH_REFLECT   = any(setFlags==24);      % reflect line
mdhBitFields.MDH_NOISEADJSCAN = any(setFlags==25);      % noise adjust scan --> Not used in NUM4
mdhBitFields.MDH_SHARENOW  = any(setFlags==26);      % all lines are acquired from the actual and previous e.g. phases
mdhBitFields.MDH_LASTMEASUREDLINE = any(setFlags==27);      % indicates that the current line is the last measured line of all succeeding e.g. phases
mdhBitFields.MDH_FIRSTSCANINSLICE = any(setFlags==28);      % indicates first scan in slice = any(setFlags==needed for time stamps)
mdhBitFields.MDH_LASTSCANINSLICE = any(setFlags==29);      % indicates  last scan in slice = any(setFlags==needed for time stamps)
mdhBitFields.MDH_TREFFECTIVEBEGIN = any(setFlags==30);      % indicates the begin time stamp for TReff = any(setFlags==triggered measurement)
mdhBitFields.MDH_TREFFECTIVEEND = any(setFlags==31);


function [magnitude,phase] = processSlices(complexSlice,acquisitionMatrixWidth)

magnitude = zeros(acquisitionMatrixWidth,acquisitionMatrixWidth,size(complexSlice,3));
phase = magnitude;
centerWindowFactor = 1/4;
    
for k = 1:size(complexSlice,3)   
    
    complexSlice(:,:,k) = echoShiftCorrection(complexSlice(:,:,k)); % center echo
    complexSlice(:,:,k) = apodization(complexSlice(:,:,k),0.5); % gibb's phenomenon, hanning window to prevent ringing
    fftdata = fftshift(ifft2(ifftshift(complexSlice(:,:,k))));
    fftdata = fftdata(:,acquisitionMatrixWidth/2+1:end-acquisitionMatrixWidth/2);
    magnitude(:,:,k) = abs(fftdata);
    
    complexSlice(:,:,k) = HighPassCmplx(complexSlice(:,:,k), centerWindowFactor);
    fftdata=fftshift(ifft2(ifftshift(complexSlice(:,:,k))));
    fftdata = fftdata(:,acquisitionMatrixWidth/2+1:end-acquisitionMatrixWidth/2);
    phase(:,:,k) = angle(fftdata);
    
end

function saveAnalyzeImage(image,filename)

temp1 = permute(image,[3 1 2]);
% height = size(image,1);
% depth = size(image,3);
analyzeWriter(temp1,1,1,1,16,filename);

%The following function is to high pass filter the phase of the cmplex data by using a window   
function filteredKspace = HighPassCmplx(cmplx_ktemp, centerWindowFactor)

[xsize,ysize] = size(cmplx_ktemp);
cmplx_ktemp3=cmplx_ktemp*0;
centerwindowsize = centerWindowFactor*min(xsize,ysize);
halfwsize=centerwindowsize/2;

cmplx_ktemp3(xsize/2-halfwsize:xsize/2+halfwsize-1, ysize/2-halfwsize:ysize/2+halfwsize-1)= ...
    cmplx_ktemp(xsize/2-halfwsize:xsize/2+halfwsize-1, ysize/2-halfwsize:ysize/2+halfwsize-1);   

cmplx_itemp2=ifft2(ifftshift(cmplx_ktemp));
cmplx_itemp22=ifft2(ifftshift(cmplx_ktemp3));

cmplx_itemp=cmplx_itemp2.*conj(cmplx_itemp22);
filteredKspace = fftshift(fft2(cmplx_itemp));

function complexSignal = echoShiftCorrection(complexSignal)

%complexSignal = ifftshift(complexSignal);
[maxVal,maxInd] = max(abs(complexSignal(:)));
[xMax,yMax] = ind2sub(size(complexSignal),maxInd);

acquisitionOneSide = size(complexSignal,1);

xshift = acquisitionOneSide/2-xMax+1;
yshift = acquisitionOneSide-yMax+1;

%Echo shift prior to any data analysis to remove global phase shift
complexSignal=circshift(complexSignal,[xshift,yshift]);
%complexSignal = fftshift(complexSignal);

function complexSignal = apodization(complexSignal,skirtingFactor)

% skirtingFactor should be between 0 and 1. The case skirtingFactor = 1
% reduces to the MATLAB hanning() function.

[xsize,ysize] = size(complexSignal);

S_h = ysize*skirtingFactor; % slope length
S_v = xsize*skirtingFactor;

skirt_h = 1/2*(1+cos((1:S_h/2)*2*pi/S_h));
skirt_v = 1/2*(1+cos((1:S_v/2)*2*pi/S_v));

f_h = [fliplr(skirt_h) ones(1,ysize-2*length(skirt_h)+1) skirt_h(1:end-1)];
f_v = [fliplr(skirt_v) ones(1,xsize-2*length(skirt_v)+1) skirt_v(1:end-1)];


complexSignal = repmat(f_v.',1,ysize).*complexSignal;
complexSignal = repmat(f_h,xsize,1).*complexSignal;
%}
