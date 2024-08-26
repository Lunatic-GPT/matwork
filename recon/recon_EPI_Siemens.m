function res = recon_EPI_Siemens(filename,varargin)
% Example:
%
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70]);
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70],'Echo',2,'Rep',50);
tic;

prefix=strtok(filename,'.');
%%{
[a,kindex]=readMeasDat(filename,'coilScaling.txt',[]);

%a=a(1:2:end-1,:)+a(2:2:end,:);
a=ifft1c(a,1);
a=a(end/4+1:end-end/4,:);
a=fft1c(a,1);

 %%
 nch=readsPar([prefix,'.pro'],'lMaximumNofRxReceiverChannels');
kindex=kindex(1:nch:end)+1;

 nro=size(a,1);
 a=reshape(a,[nro,nch,size(a,2)/nch]);
 a=permute(a,[1,3,2]);
 npe=readsPar([prefix,'.pro'],'lPhaseEncodingLines');
 nrefpe=readsPar([prefix,'.pro'],'lRefLinesPE');
accel=readsPar([prefix,'.pro'],'lAccelFactPE');
nrep=readsPar([prefix,'.pro'],'lRepetitions');
if isempty(nrep)  % parameter lRepetitions not store in .pro file when nrep =0
    nrep=0;
end
 nrep=nrep+1;
%%
d_nois=a(:,1,:);
d_cal=a(:,2:1+(nrefpe/accel+3)*accel,:);
d_cal=reshape(d_cal,[nro,(nrefpe/accel+3),accel,nch]);
%d_cal(:,1,:,:)=mean(d_cal(:,[1,3],:,:),2);
if mod(size(d_cal,2),2)==0
 d_cal(:,3,:,:)=[];
else
 d_cal(:,1,:,:)=[];
end

d_cal(:,1:2:end,:,:)=flipdim(d_cal(:,1:2:end,:,:),1);
d_cal_pc=epi_phase_corr(d_cal(:,3:end,:,:),(d_cal(:,[1,2],:,:)));

%%
%kindex=load([prefix,'.dat_kindex.log']);
k_cal=kindex(2:1+(nrefpe/accel+3)*accel);
k_cal=reshape(k_cal,[nrefpe/accel+3,accel]);
k_cal=k_cal(4:end,:);
k_cal=k_cal(:);
k_cal=k_cal-min(k_cal)+1;

d_cal_pc2=reshape(d_cal_pc,[nro,nrefpe,nch]);

d_cal_pc2(:,k_cal,:)=d_cal_pc2;

%%
k_data=kindex(2+(nrefpe/accel+3)*accel:end);
k_data=reshape(k_data,[length(k_data)/nrep,nrep]);
k_data=k_data(4:end,1);
k_center=kindex(2);
nk_data=length(k_data);
d_data=a(:,2+(nrefpe/accel+3)*accel:end,:);
d_data=reshape(d_data,[nro,nk_data+3,nrep,nch]);

%d_data(:,1,:,:)=mean(d_data(:,[1,3],:,:),2);
if mod(size(d_data,2),2)==0
 d_data(:,3,:,:)=[];
else
 d_data(:,1,:,:)=[];
end

if mod(size(d_data,2),2)~=mod(size(d_cal,2),2)
   d_cal_pc2=flip(d_cal_pc2,1); 
end
%d_data(:,1,:,:)=[];
d_data(:,1:2:end,:,:)=flipdim(d_data(:,1:2:end,:,:),1);

d_data_pc=epi_phase_corr(d_data(:,3:end,:,:),(d_data(:,[1,2],:,:)));

d_data_pc2=zeros(nro,max(k_data)+accel-1,nrep,nch);
d_data_pc2(:,k_data,:,:)=d_data_pc;

%}

Rkall=recon_grappa(d_data_pc2,d_cal_pc2,k_data);
Rkall=Rkall(:,1:end-2,:,:);
[img,ph]=partialFT(Rkall,k_center-1);
a2=sqrt(sum(abs(img).^2,4));
img=reshape(a2,[nro,npe,1,nrep]);


save(['recon_',prefix],'img','ph');

fprintf('Read data took %s s\n',toc);


function Rkall=recon_grappa(data,dref,isamp)


 nro=size(data,1);
 npe=size(data,2);
 nrep=size(data,3);
 nch=size(data,4);
 
 Rkall=zeros(nro,npe,nrep,nch);
  
  kernelX=[-2,-1,0,1,2];
 kernY=[-4,-1,2,5];
    [tmp,coef1]=myGRAPPA_xp(dref,1:size(dref,2),1:size(dref,2),kernY,kernelX,0);
    
    
     kernY=[-5,-2,1,4];
     [tmp,coef2]=myGRAPPA_xp(dref,1:size(dref,2),1:size(dref,2),kernY,kernelX,0);
    
     
    for irep=1:nrep
    disp(irep);
    %    data=fft1c(squeeze(dref3(:,:,isl,:)),1);
   
    kernY=[-4,-1,2,5];
    
    Rk=myGRAPPA_xp(squeeze(data(:,:,irep,:)),isamp,0,kernY,kernelX,coef1);
    
     kernY=[-5,-2,1,4];
    
     Rk2=myGRAPPA_xp(squeeze(data(:,:,irep,:)),isamp,0,kernY,kernelX,coef2);
    
     Rk(:,isamp+2,:)=Rk2(:,isamp+2,:);
     Rkall(:,:,irep,:)=Rk;
    end
    
 %   img1b=ifft2c(Rkall);
      
%img2=sqrt(sum(abs(img1b).^2,4));

   

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

