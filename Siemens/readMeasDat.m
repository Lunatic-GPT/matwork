function [sdata,snavData,dummyData] = readMeasDat(fname,max_nlines,offset,remove_os)
% fname: data file name
% max_nlines; the max number of lines to read
% offset: the starting offset; 0 start from the first line.
% remove_os: If true, remove the 2* oversampling along readout.  If
% remove_os is true, then actualData is already Inv Fourier
% transformed along readout. default: false.
% 
% sdata and snavData has the following fields:
% Line, Partition, Set, Slice, Channel, Flags, Repetition, Data

% Line and Partition are 0 based;
%
%  10/31/2019: format changed; output flags
% [actualData,ushLine,ushPartition,ushSlice,ushSet,PMUTimeStamp,icePara,navData,ushLineNav,ushPartitionNav,ushSliceNav,RepetitionNav,dummyData] = readMeasDat(fname,max_nlines,offset,remove_os)
% 

tic;

if ~exist('remove_os','var')
    remove_os=false;
end

% VB17 RaidFile structure

fprintf('Reading file %s\n',fname);
% FILE *in = fopen("D:\\temp\\meas_orig.dat", "rb");
fid = fopen(fname,'rb');
% int32_t nProtHeaderLen = 0; // read the protocol header length
% fread(&nProtHeaderLen, 1, 4, in);
nProtHeaderLen = double(fread(fid,1,'*int32'));
% fseek(in, nProtHeaderLen - 4, SEEK_CUR);

%nProtHeaderLen=10000;
[tmp,xprotocol]=fileparts(strtok2(fname,'.'));
xprotocol=[xprotocol,'.xprotocol'];

prefix=strtok2(fname,'.');

if ~exist(prefix,'dir')
    mkdir(prefix);
end

if exist(fullfile(prefix,xprotocol),'file')
    fseek(fid, nProtHeaderLen - 4, 0);
else
    header=fread(fid,nProtHeaderLen-4,'*char');
    fid2=fopen(fullfile(prefix,xprotocol),'w');
    
    for i=5:length(header)-12
        fprintf(fid2,'%c',header(i));
    end
    fclose(fid2);
    
end

if ~exist(fullfile(prefix,[filename(prefix),'.pro']),'file')
    extp(fullfile(prefix,xprotocol));
end



if ~exist('offset','var')
    offset=0;
end
if ~exist('max_nlines','var')
    max_nlines=1000;
end


acqEnded = 0;
mdhCounter = 1;

% fidLog = fopen('log.txt','wt');
oldMostSig=-1;

Data=single([]);
navData=single([]);
dummyData=single([]);
DataLine=uint16([]);
DataPartition=uint16([]);
DataSlice=uint16([]);
DataCha=uint16([]);
DataRepetition=uint16([]);
DataSet=uint16([]);
DataFlag=uint32([]);


navLine=uint16([]);
navPartition=uint16([]);
navSlice=uint16([]);
navCha=uint16([]);
navFlag=uint32([]);
navRepetition=uint16([]);
navSet=uint16([]);


icePara=uint16([]);
freePara=uint16([]);
iceParaNav=uint16([]);
freeParaNav=uint16([]);

PMUTimeStamp=uint32([]);
PMUTimeStampNav=uint32([]);

TimeStamp=uint32([]);
TimeStampNav=uint32([]);

RepetitionNav=uint16([]);
 
while 1%~acqEnded
    
    
    % MDH Header Information Retrieval Starts Here...
    ulFlagsAndDMALength = fread(fid,1,'*uint32');
    if isempty(ulFlagsAndDMALength)
       disp('isempty(ulFlagsAndDMALength)');
        break;
    end
    lMeasUID = fread(fid,1,'*int32');
    ulScanCounter = fread(fid,1,'*uint32');
    tmpTimeStamp = fread(fid,1,'*uint32');
    tmpPMUTimeStamp = fread(fid,1,'*uint32');
    
    aulEvalInfoMaskMostSig = fread(fid,1,'*uint32');%evaluation info mask field,  first part
    aulEvalInfoMaskLeastSig = fread(fid,1,'*uint32');%evaluation info mask field, second part
    mdhBitFields = determineBitFields(aulEvalInfoMaskMostSig);
    
    
    
    if oldMostSig ~=aulEvalInfoMaskMostSig
        oldMostSig=aulEvalInfoMaskMostSig;
        %   mdhBitFields = determineBitFields(aulEvalInfoMaskMostSig);
        acqEnded = mdhBitFields.MDH_ACQEND;
        doscale = mdhBitFields.MDH_RAWDATACORRECTION;
        
    
    end
    
 
    % acqEnd is automatically added (i.e if not by the sequence); In this case, if not
    % skipped here (by "if acqEnded") will generate an exception below because the column length
    % mismatches (numberOfColumns=16);
    
    if feof(fid) %||acqEnded    
        PMUTimeStamp(end)=[];
        TimeStamp(end)=[];
        fprintf('feof(fid) is true'); 
        break;
    end
    ushSamplesInScan = fread(fid,1,'*uint16'); % # of samples acquired in scan
    ushUsedChannels = fread(fid,1,'uint16');
    ushLine_tmp = fread(fid,1,'*uint16');                  % line index                   */
    ushAcquisition = fread(fid,1,'*uint16');           % acquisition index            */
    ushSlice_tmp = fread(fid,1,'*uint16');                 % slice index                  */
    ushPartition_tmp = fread(fid,1,'*uint16');             % partition index              */
    ushEcho = fread(fid,1,'*uint16');                  % echo index                   */
    ushPhase = fread(fid,1,'*uint16');                 % phase index                  */
    Repetition_tmp = fread(fid,1,'*uint16');            % measurement repeat index     */
    ushSet=fread(fid,1,'*uint16');
    fseek(fid,6*2,0); % Other sLoopCounter variables
    fseek(fid,4,0); % sCutOff variables
    ushKSpaceCenterColumn = fread(fid,1,'*uint16');
    ushCoilSelect = fread(fid,1,'*uint16');
    fReadOutOffCenter = fread(fid,1,'*float32');
    ulTimeSinceLastRF = fread(fid,1,'*uint32');
    ushKSpaceCenterLineNo = fread(fid,1,'*uint16');
    ushKSpaceCenterPartitionNo = fread(fid,1,'*uint16');
    %   fseek(fid,8,0); % aushIceProgramPara
    tmpIcePara=fread(fid,4,'*uint16');
    %fseek(fid,8,0); % aushicePara
    tmpFreePara=fread(fid,4,'*uint16');

    %fseek(fid,28,0); % sSliceData -> replaced with fSag ... aflQuaternion
    fSag = fread(fid,1,'*float32');
    fCor = fread(fid,1,'*float32');
    fTra = fread(fid,1,'*float32');
    aflQuaternion = fread(fid,4,'*float32');
    %disp([fSag fCor fTra aflQuaternion'])
    ushChannelID = fread(fid,1,'*uint16');
    fseek(fid,2,0);
    % MDH Header Information Retrieval Ends Here...
    numberOfColumns = double(ushSamplesInScan);
      
    if (mdhCounter==1) &&offset>0
        fseek(fid,2*numberOfColumns*4*offset+(offset-1)*128,0);
        offset=0;
        continue;
    end
    %  if exist('chan','var') &&~isempty(chan)&& (ushChannelID+1)~=chan
    %        fseek(fid,2*ushSamplesInScan*4,0);
    %  else
    try
        tmpData = fread(fid,2*numberOfColumns,'*float32');
    catch
        
        disp('ok');
    end
    
        rl=tmpData(1:2:end);
        im=tmpData(2:2:end);
        
    if doscale
        
        if ~exist('scl','var')
            if ushUsedChannels>1
              try
                scl=readCoilScalingFromProtocol(fullfile(prefix,xprotocol),ushUsedChannels);
                tmpData=(rl+1i*im)*scl(ushChannelID+1);
             catch
                fprintf('error reading scale; number of used channels = %d; mdhcounter = %d; LineNo = %d; ParNo = %d\n',...
                      ushUsedChannels,mdhCounter,ushKSpaceCenterLineNo,ushKSpaceCenterPartitionNo);
                tmpData=(rl+1i*im);

                end
            else
                tmpData=(rl+1i*im);
            end
        end
        
    else
        
        tmpData=(rl+1i*im);
    end
    
    
    if (mdhBitFields.MDH_REFPHASESTABSCAN || mdhBitFields.MDH_RTFEEDBACK)
        try
         
        %  fprintf('Nav data found.');
          count=size(navData,2);
      
         if ~remove_os
            navData(:,count+1)=tmpData;
         else
            ftmpData=ifft1c(tmpData,1);
            len=length(ftmpData);
            navData(:,count+1)=ftmpData(len/4+1:end-len/4);    
         end
         navRepetition(count+1)=Repetition_tmp;
         navLine(count+1)=ushLine_tmp;
         navPartition(count+1)=ushPartition_tmp;
         navSlice(count+1)=ushSlice_tmp;
         navCha(count+1)=ushChannelID;
         navSet(count+1)=ushSet;
         navFlag(count+1)=aulEvalInfoMaskMostSig;
         TimeStampNav(count+1) = tmpTimeStamp;
         PMUTimeStampNav(count+1) = tmpPMUTimeStamp;
         iceParaNav(count+1,:)=tmpIcePara;
         freeParaNav(count+1,:)=tmpFreePara;
          
        catch 
            disp('error reading nav data');
        end
    elseif (mdhBitFields.MDH_PHASESTABSCAN)
         
        count=size(dummyData,2);    
        if ~remove_os
            dummyData(:,count+1)=tmpData;
        else
            ftmpData=ifft1c(tmpData,1);
            len=length(ftmpData);
            dummyData(:,count+1)=ftmpData(len/4+1:end-len/4);    
        end
        
    else
       
             count=size(Data,2);
        
        try
            if ~remove_os
                Data(:,count+1)=tmpData;
            else
                ftmpData=ifft1c(tmpData,1);               
                len=length(ftmpData);
                Data(:,count+1)=ftmpData(len/4+1:end-len/4);
            end
        %if count<68000
           % fprintf('setting data: size Data = %d; tmpData = %d, ch = %d, line = %d, part = %d\n',size(Data,1),size(tmpData,1),ushChannelID,ushLine_tmp,ushPartition_tmp);
        %end
        catch
           if ushChannelID==0
             fprintf('error setting data: size Data = %d; tmpData = %d, ch = %d, line = %d, part = %d\n',size(Data,1),size(tmpData,1),ushChannelID,ushLine_tmp,ushPartition_tmp);
           end
            %break;
  
        end
            
        DataLine(count+1)=ushLine_tmp;
        DataPartition(count+1)=ushPartition_tmp;
        DataSlice(count+1)=ushSlice_tmp;
        DataSet(count+1)=ushSet;
        DataCha(count+1)=ushChannelID;
        DataFlag(count+1)=aulEvalInfoMaskMostSig;
        DataRepetition(count+1)=Repetition_tmp;
        TimeStamp(count+1) = tmpTimeStamp;
        PMUTimeStamp(count+1) = tmpPMUTimeStamp;
        icePara(count+1,:)=tmpIcePara;
        freePara(count+1,:)=tmpFreePara;
    
    end
        
    
    
    if mdhCounter==max_nlines
       fprintf('max line reached');
        break;
    end
    
    mdhCounter = mdhCounter+1;
    
    if mod(mdhCounter,1000)==0
        disp(mdhCounter);
    end
end
fclose(fid);


sdata.Line=DataLine;
sdata.Partition=DataPartition;
sdata.Set=DataSet;
sdata.Slice=DataSlice;
sdata.Channel=DataCha;
sdata.Flags=DataFlag;
sdata.Repetition=DataRepetition;
sdata.Data=Data;
sdata.PMUTimeStamp = PMUTimeStamp;
sdata.TimeStamp = TimeStamp;
sdata.icePara=icePara;
sdata.freePara=freePara;
    

snavData.Data=navData;
snavData.Line=navLine;
snavData.Partition=navPartition;
snavData.Set=navSet;
snavData.Slice=navSlice;
snavData.Channel=navCha;
snavData.Flags=navFlag;
snavData.Repetition=navRepetition;
snavData.PMUTimeStamp = PMUTimeStampNav;
snavData.TimeStamp = TimeStampNav;
snavData.icePara=iceParaNav;
snavData.freePara=freeParaNav;


% ushLine=ushLine(1:size(actualData,2));
% ushPartition=ushPartition(1:size(actualData,2));
% ushSlice=ushSlice(1:size(actualData,2));
fprintf('Read data took %s s\n',toc);

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

function saveAnalyzeImage(image,fname)

temp1 = permute(image,[3 1 2]);
% height = size(image,1);
% depth = size(image,3);
analyzeWriter(temp1,1,1,1,16,fname);

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

