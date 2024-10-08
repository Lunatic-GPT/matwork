function [kSpaceData,noiseSamples,sampledLineIndices,accelerationFactor] = readMeasDat_NavCorr3D(filename,mtrx,channel)
% Example:
%
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70]);
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70],'Echo',2,'Rep',50);


% VB17 RaidFile structure

% FILE *in = fopen("D:\\temp\\meas_orig.dat", "rb");
fid = fopen(filename,'rb');
% int32_t nProtHeaderLen = 0; // read the protocol header length
% fread(&nProtHeaderLen, 1, 4, in);
nProtHeaderLen = double(fread(fid,1,'*int32'));

if exist([filename,'_protocol'],'file')
   fseek(fid, nProtHeaderLen - 4, 0);
else
  header=fread(fid,nProtHeaderLen-4,'*char');
  fid2=fopen([filename,'_protocol'],'w');

  for i=5:length(header)-12
   fprintf(fid2,'%c',header(i));
  end
  fclose(fid2);

end
%fseek(fid,nProtHeaderLen-4,0);

amplitudeScaleFactor = 80 * 20 * 131072 / 65536;

amplitudeScaleFactor = amplitudeScaleFactor*20000;

acqEnded = 0;
mdhCounter = 1;
sliceIndices = [];
sampledLineIndices = [];

noiseSamples = zeros(32,256);

while ~acqEnded
   % MDH Header Information Retrieval Starts Here...
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
    fseek(fid,28,0); % sSliceData
    ushChannelID = fread(fid,1,'*uint16');
    fseek(fid,2,0);
    % MDH Header Information Retrieval Ends Here...   
    mdhBitFields = determineBitFields(aulEvalInfoMaskMostSig);   
    acqEnded = mdhBitFields.MDH_ACQEND;
    actualData = fread(fid,2*(double(ushSamplesInScan)),'*float32');
    realPart = actualData(1:2:end-1);%*amplitudeScaleFactor;
    imaginaryPart = actualData(2:2:end);%*amplitudeScaleFactor;
    complexData = realPart+1i*imaginaryPart;
    
    
    
    if mdhBitFields.MDH_NOISEADJSCAN
        if sum(abs(noiseSamples(:)))==0
            noiseSamples = zeros(32,length(complexData));
        end
        noiseSamples(ushChannelID+1,:) = complexData;
        continue;
    end    
    
    if acqEnded
        %return;
        break;
    end
    
    if mdhCounter==1 && exist('mtrx','var') && ~isempty(mtrx)
        kSpaceData=single(zeros([mtrx,length(channel)]));
    end
    if mdhBitFields.MDH_REFLECT
        complexData = flipud(complexData); 
    end
    
    if any(ushChannelID+1==channel)
        
     % kSpaceData(ushLine+1,:,ushPartition+1,ushChannelID+1,ushEcho+1) = complexData;
      ichan=find(ushChannelID+1==channel);
      kSpaceData(:,ushLine+1,ushPartition+1,ichan) = complexData;
      
    end
    
    currentLine = ushLine+1;
    currentAcquisition = ushAcquisition+1;
    currentSlice = ushSlice+1;
    currentRepetition = ushRepetition+1;
    currentChannel = ushChannelID+1; 
    currentEcho = ushEcho+1;
%     
    if ~any(sampledLineIndices==currentLine) && ~mdhBitFields.MDH_NOISEADJSCAN 
        sampledLineIndices(end+1)=currentLine;        
    end
    
    if currentChannel==1 && ushPartition==0 && ushEcho == 0       
        disp(['Channel: ' num2str(currentChannel)...
            ' Slice: ' num2str(currentSlice)...
            ' Line: ' num2str(currentLine) ...
            ' Partition: ' num2str(ushPartition+1) ...
            ' Echo: ' num2str(currentEcho) ...
            ' # Samples: ' num2str(ushSamplesInScan)]);                
      if currentLine ==190 
          disp('pause');
      end
    end
    
    mdhCounter = mdhCounter+1;   
    
end
fclose(fid);
sampledLineIndices = sort(sampledLineIndices);
accelerationFactor = diff(sampledLineIndices(1:2));


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

