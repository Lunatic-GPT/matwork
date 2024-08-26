function res = readMeasDat_mbDiff(filename,varargin)
% Example:
%
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70]);
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70],'Echo',2,'Rep',50);
tic;
numberOfPE = 0;
numberOfSlices = 0;
numberOfEchoes = 1;
numberOfRepetitions = 1;

if nargin >= 2
    numberOfPE = varargin{1}(1);    
    numberOfSlices = varargin{1}(3);   
    if nargin == 4 && strcmpi('echo',varargin{2})
        numberOfEchoes = varargin{3};
    end
    if nargin == 6 && strcmpi('rep',varargin{4})
        numberOfRepetitions = varargin{5};
    end
end

% VB17 RaidFile structure

% FILE *in = fopen("D:\\temp\\meas_orig.dat", "rb");
fid = fopen(filename,'rb');
% int32_t nProtHeaderLen = 0; // read the protocol header length
% fread(&nProtHeaderLen, 1, 4, in);
nProtHeaderLen = double(fread(fid,1,'*int32'));
% fseek(in, nProtHeaderLen - 4, SEEK_CUR);

if exist([filename,'_protocol'],'file')
   fseek(fid, nProtHeaderLen - 4, 0);
else
  header=fread(fid,nProtHeaderLen-4,'*char');
  fid2=fopen([filename,'_protocol'],'w');

  for i=5:length(header)-12
   fprintf(fid2,'%c',header(i));
  end
  fclose(fid2);

  extract_protocol([filename,'_protocol']);
end


amplitudeScaleFactor = 1;

acqEnded = 0;
mdhCounter = 1;
sliceIndices = [];
sampledLineIndices = [];

% fidLog = fopen('log.txt','wt');
oldMostSig=-1;
    
    ulFlagsAndDMALength = fread(fid,1,'*uint32');
    lMeasUID = fread(fid,1,'*int32');
    ulScanCounter = fread(fid,1,'*uint32');
    ulTimeStamp = fread(fid,1,'*uint32');
    ulPMUTimeStamp = fread(fid,1,'*uint32');
    aulEvalInfoMaskMostSig = fread(fid,1,'*uint32');%evaluation info mask field,  first part    
    aulEvalInfoMaskLeastSig = fread(fid,1,'*uint32');%evaluation info mask field, second part    
    ushSamplesInScan = fread(fid,1,'*uint16'); % # of samples acquired in scan 
    ushUsedChannels = fread(fid,1,'uint16');
    ushLine = fread(fid,1,'*uint16');                  % line index                   */
    ushAcquisition = fread(fid,1,'*uint16');           % acquisition index            */
    ushSlice = fread(fid,1,'*uint16');                 % slice index                  */
    ushPartition = fread(fid,1,'*uint16');             % partition index              */
    ushEcho = fread(fid,1,'*uint16');                  % echo index                   */	
    ushPhase = fread(fid,1,'*uint16');                 % phase index                  */
    ushRepetition = fread(fid,1,'*uint16');            % measurement repeat index     */    
    fseek(fid,7*2,0); % Other sLoopCounter variables
    fseek(fid,4,0); % sCutOff variables
    ushKSpaceCenterColumn = fread(fid,1,'*uint16');
    ushCoilSelect = fread(fid,1,'*uint16');
    fReadOutOffCenter = fread(fid,1,'*float32');
    ulTimeSinceLastRF = fread(fid,1,'*uint32');
    ushKSpaceCenterLineNo = fread(fid,1,'*uint16');
    ushKSpaceCenterPartitionNo = fread(fid,1,'*uint16');
    fseek(fid,8,0); % aushIceProgramPara
    fseek(fid,8,0); % aushFreePara
    %fseek(fid,28,0); % sSliceData -> replaced with fSag ... aflQuaternion
    fSag = fread(fid,1,'*float32');
    fCor = fread(fid,1,'*float32');
    fTra = fread(fid,1,'*float32');
    aflQuaternion = fread(fid,4,'*float32');    
    %disp([fSag fCor fTra aflQuaternion'])
    ushChannelID = fread(fid,1,'*uint16');
    fseek(fid,2,0);
    
    fseek(fid,-128,0);
    % MDH Header Information Retrieval Ends Here...
    
    nro=double(ushSamplesInScan);
    prefix=strtok(filename,'.');
    npe=readsPar([prefix,'.pro'],'lPhaseEncodingLines');
    
    nb=readsPar([prefix,'.pro'],'lBCSeqExcitationMode');
    % Slices    

 
    nch=double(ushUsedChannels);

    nsl=readsPar([prefix,'.pro'],'lFinalMatrixSizeSlice');
    nref=readsPar([prefix,'.pro'],'lRefLinesPE');
    accel=readsPar([prefix,'.pro'],'lAccelFactPE');
    ndiff=readsPar([prefix,'.pro'],'lDiffDirections');
    
   
    
    lnois=nch*(nro*2+32);  %k=0
   
    lref=(nro*2+32)*nch*(nref+3*accel)*nsl;
    l1b=(nro*2+32)*(npe/accel-1+3)*nsl*nch;
    lmb=(nro*2+32)*(npe/accel-1+3)*nsl/nb*nch;
    
    res=fread(fid,lnois,'*float32');
    
    dnois=reshape(res,[(nro*2+32),nch]);
    dnois=dnois(33:end,:);
    dnois=dnois(1:2:end,:)+1i*dnois(2:2:end,:);
    
    
   
   
    res=fread(fid,lref,'*float32');
    dref=reshape(res,[nro*2+32,nch,nref/accel+3,nsl,accel]);
    dref=dref(33:end,:,:,:,:);
    dref=dref(1:2:end,:,:,:,:)+1i*dref(2:2:end,:,:,:,:);
    dref2=permute(dref(:,:,1:end,:,:),[1,5,3,4,2]);
  
    if mod(nsl,2)==0
     sorder=[2:2:nsl,1:2:nsl-1];
    else
     sorder=[1:2:nsl-1,2:2:nsl];
    end
    dref2(:,:,:,sorder,:)=dref2;
    
    %dreftmp=dref2;
   % dref2(:,:,1,:,:)=dref2(:,:,3,:,:); 
    dref2(:,:,1,:,:)=[];
    dref2(:,:,1:2,:,:)=dref2(:,:,2:-1:1,:,:);
   
    
    dref2(:,:,2:2:end,:,:)=dref2(end:-1:1,:,2:2:end,:,:);
  
   dref2=ifft1c(dref2,1);
     
   ph=angle(dref2(:,:,1:2,:,:));
   dref2=dref2(:,:,3:end,:,:);
   
   dref2(:,:,1:2:end,:,:)=dref2(:,:,1:2:end,:,:).*repmat(exp(-1i*ph(:,:,1,:,:)),[1,1,size(dref2,3)/2,1,1]);
   dref2(:,:,2:2:end,:,:)=dref2(:,:,2:2:end,:,:).*repmat(exp(-1i*ph(:,:,2,:,:)),[1,1,size(dref2,3)/2,1,1]);
   
    dref3=reshape(dref2,[nro,nref,nsl,nch]);
    
    kernelY=-3:3;
    kernelX=-3:3;
    fdref3=fft1c(dref3,1);
 %   [Sn,coef] = myGRAPPA_Sense_xp(squeeze(fdref3(:,:,44,:)),1:size(fdref3,2),1:size(fdref3,2),kernelY,kernelX,false)

%    img=ifft1c(dref3,2);
%   figure;
% img2=sqrt(sum(abs(img).^2,4));
% 
% for i=1:52
%     
%     subplot(7,8,i);
%     shft=floor((i-1)/13);
%     tmp=circshift(img2(65:192,:,i),[0,-shft*21]);
%     imshow(fliplr(flipud(tmp)),[]);
% end

%%
    res=fread(fid,l1b,'*float32');
   
    % fclose(res);
       
    d1b=reshape(res,[(nro*2+32),nch,npe/accel+3-1,nsl]);
    d1b=d1b(33:end,:,:,:);
    d1b=d1b(1:2:end,:,:,:)+1i*d1b(2:2:end,:,:,:);
    
    d1b=permute(d1b,[1,3,4,2]);
    if mod(nsl,2)==0
     sorder=[2:2:nsl,1:2:nsl-1];
    else
     sorder=[1:2:nsl-1,2:2:nsl];
    end
    d1b(:,:,sorder,:)=d1b;
   
    d1b(:,1,:,:)=[];
    d1b(:,1:2,:,:)=d1b(:,2:-1:1,:,:);
    
    d1b(:,2:2:end,:,:)=d1b(end:-1:1,2:2:end,:,:);
    
    ref1b=d1b(:,1:2,:,:);
    d1b=d1b(:,3:end,:,:);
  
%% uncomment to recon single band images.    
    d1bc=phase_corr(d1b,ref1b);
    d1b2=zeros(nro,npe,nsl,nch);

    isamp=2:accel:npe-2;
       
     d1b2(:,isamp,:,:)=d1bc;   
    
    % d1bk=flipdim(d1b2,1);
   
   recon_grappa(d1b2(:,:,44,:),dref3(:,:,44,:),isamp,1);
  %%
 
   
   
    method='slice_grappa';
switch method 
  case 'slice_grappa'
         
    kernY=-3:3;
    kernX=-3:3;
             
 
     clear coef;
     
    for j=1:nb
        for i=1:nsl/nb 
            disp(i+(j-1)*nsl/nb);
         coef(:,:,i+(j-1)*nsl/nb)=myGRAPPA_2data_xp(squeeze(sum(d1b(:,:,i:end/nb:end,:),3)),squeeze(d1b(:,:,i+(j-1)*nsl/nb,:)),1:size(d1b,2),1:size(d1b,2),kernY,kernX,0);  
        end
    end
    save('coef','coef');
    
   for idiff=1:ndiff+1
     res=fread(fid,lmb,'*float32');
   
       
    dmb=reshape(res,[(nro*2+32),nch,npe/accel+3-1,nsl/nb]);
    dmb=dmb(33:end,:,:,:,:);
    dmb=dmb(1:2:end,:,:,:,:)+1i*dmb(2:2:end,:,:,:,:);
    
    dmb=permute(dmb,[1,3,4,2,5]);
    
    
     
     
    sorder=(0:nsl/nb-1)*5;
    sorder=mod(sorder,22)+1;
    dmb=dmb(:,:,sorder,:);
   
    dmb(:,1,:,:,:)=[];
    
    dmb(:,1:2,:,:,:)=dmb(:,2:-1:1,:,:,:);
    dmb(:,2:2:end,:,:,:)=dmb(end:-1:1,2:2:end,:,:,:);
    
    refmb=dmb(:,1:2,:,:,:);
    
    dmb=dmb(:,3:end,:,:,:);
    
     dmb4=zeros(nro,npe,nsl,nch);
              isamp=2:2:126;  
    for j=1:nb
        for i=1:nsl/nb
         dmb_sep=myGRAPPA_2data_xp(squeeze(dmb(:,:,i,:,1)),[],1:size(dmb,2),1:size(dmb,2),kernY,kernX,coef(:,:,i+(j-1)*nsl/nb));
         dmb_sep=reshape(dmb_sep,[size(dmb_sep,1),size(dmb_sep,2),1,size(dmb_sep,3)]); 
         dmb_sep2=phase_corr(dmb_sep,ref1b(:,:,i+(j-1)*nsl/nb,:));

         dmb4(:,isamp,i+(j-1)*nsl/nb,:)=dmb_sep2;%flipdim(dmb_sep2,1);   
        end
    end
   
     recon_grappa(dmb4,dref3,isamp,idiff);
   end
   
    fclose(fid);
     
  case 'cat_slice'
        
    for ir=1
        
    dmb2=squeeze(dmb(:,:,:,:,ir));
                                                  
    fd1b=fft(d1b,[],3);
    fd1b=circshift(fd1b,[0,0,ceil(nsl/2),0]);
    
    dmb3=repmat(zeros(size(dmb2)),[1,1,4,1]);
    ismp=1:4:size(dmb3,3);
    dmb3(:,:,1:4:end,:)=fft(dmb2,[],3);
    
    kernY=[-5,-1,3,7];
    kernX=-2:2;
    
    len=size(d1b,2);
    
    %[tmp,tmp2,coef]=myGRAPPA_xp(squeeze(fd1b(:,round(end/2),:,:)),1:nsl,1:nsl,kernY,kernX,0);
    
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

     isamp=2:2:126;
        
      dmb4(:,isamp,:,:)=dmb_sep;   
     
   %   dmb4=flipdim(dmb4,1);
      recon_grappa(dmb4(:,:,44:45,:),dref3(:,:,44:45,:),isamp);
 
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

function recon_grappa(d1bk,dref3,isamp,idiff)


 nro=size(d1bk,1);
 npe=size(d1bk,2);
 nsl=size(d1bk,3);
 nch=size(d1bk,4);
 
 Rkall=zeros(nro,npe,nsl,nch);
  
    for isl=1:nsl
    disp(isl);
        data=fft1c(squeeze(dref3(:,:,isl,:)),1);
   
    kernY=[-3,-1,1,3];
    [tmp,coef]=myGRAPPA_xp(data,1:size(dref3,2),3:size(dref3,2)-2,kernY,[-1,0,1],0);
    
    
    [tempRk]=myGRAPPA_xp(squeeze(d1bk(:,:,isl,:)),isamp,0,kernY,[-1,0,1],coef);
    
     Rkall(:,:,isl,:)=Rk;
    end
    
   


%figure;
for i=1:nsl
    
 % subplot(8,11,i);
    shft=floor((i-1)/22);
  shft=1;
    Rkall2(:,:,i,:)=fov_shift_kspace(Rkall(:,:,i,:),[0,(shft)*43],[size(Rkall,1),size(Rkall,2)]);
    
    %img2(:,:,i)=circshift(img2(:,:,i),[0,(-shft)*43]);
  %  tmp=abs(img2(:,:,i));
  %  imshow(fliplr((tmp)),[]);drawnow;
end 

 img1b=ifft2c(Rkall2);
 img2=sqrt(sum(abs(img1b).^2,4));
img2=img2(65:192,:,:);

Rkall=single(Rkall2);
img2=single(img2);
 %  save(sprintf('recon_grappa_diff256dir_%03d',idiff),'img2','Rkall');
  save(sprintf('recon_grappa_singleband_256dir_%03d',idiff),'img2','Rkall');

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

