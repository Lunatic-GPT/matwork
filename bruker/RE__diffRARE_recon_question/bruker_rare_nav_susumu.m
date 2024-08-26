% this program is only for signle channel receiver data use
% there is another program "bruker_rare_nav_susumu_multiReceiver.m" for
% multireceiver

% The purpose of this function is to seperate the rare Raw data and 
% two navigator echo data into two seperate files. The rare raw data are 
% stored into file named "fid"; The navigator echo data are stored into 
% "fid.nav"
% The second function is to  do correctio according to Dr. Susumu Mori's paper
% MRM 40:511-516(1998) before writing rare raw data 
function bruker_rare_nav_susumu(fileName)

% modify following parameters for your scan
rareFactor = 4;
navFactor = 2;  % means twin echoes.
readSize = 205;
phaseSize = 104;
phase2Size = 8;  
diffSize = 4;   % # of b0 image plus # of diffusion directions. 
if navFactor == 0
    flagSusumu = 0;
else
    flagSusumu = 1;
end;
flagSusumu = 1;  %   turn on/off correction


%******* the following code needs not to be modified ******************
dataFormat = 'ieee-le';
flagKFormat = true; 
flagPlotPhase = 0;
 
 
fid = fopen(fileName, 'r');
if fid == -1
    error( sprintf( 'file %s can not open ', fid));
    return;
end;

fileNameRare = strcat(fileName, '_rare');
fileNameNav = strcat(fileName, '_nav');

fid_rare = fopen(fileNameRare, 'w');
if fid_rare == -1
    error( sprintf( 'file %s can not open ', fid_rare));
    return;
end;

fid_nav = fopen(fileNameNav, 'w');
if fid_nav == -1
    error( sprintf( 'file %s can not open ', fid_nav));
    return;
end;

position_fid = ftell(fid);
position_rare = ftell(fid_rare);
position_nav = ftell(fid_nav);

% aa = fread(fid, 'int32');

if flagKFormat == true; 
    blockSize = ceil((readSize*2*4)/1024);
    skipSize = blockSize*(1024/4)-readSize*2 ;
    % readSize = blockSize*(1024/4)/2;
end;

%  tempTest =  fread(fid, readSize*2+skipSize, 'int32', 'ieee-le'); 
% tempTest1 = tempTest(1:2:end) + j*tempTest(2:2:end);
dummy_counter(1:diffSize) = 1;
segmentSize = (rareFactor+navFactor)*phaseSize/rareFactor;


for kk = 1:phase2Size
    for jj = 1:(rareFactor+navFactor):(phaseSize*(rareFactor+navFactor)/rareFactor)
        for llDiff = 1:diffSize
            posDataRead = 0;
            for ii = 1:(rareFactor+navFactor)
                if ii>rareFactor         
                   tempBufferNav(1:readSize*2) = fread(fid, readSize*2, 'int32', dataFormat); 
                  if (mod(ii,2) == 1)  
                      phaseOdd = tempBufferNav(1:2:end) + 1i*tempBufferNav(2:2:end);
                  else
                      phaseEven = tempBufferNav(1:2:end) + 1i*tempBufferNav(2:2:end);
                  end;
                  
                 count =  fwrite(fid_nav, tempBufferNav, 'int32', dataFormat);
                   % display(sprintf('the current file navigator echo position is %d  \n', position_nav));
                    if count ~= readSize*2
                        error(' writing error into navigator fid file ');
                        return;
                    end;
                    if flagKFormat == true 
                       grabage = fread(fid, skipSize , 'int32', dataFormat); 
                    end;
                   % fid_position = ftell(fid);
                   % display(sprintf('the current file position after skip garbage is %d  \n', fid_position));
                    % fid_position = fid_position + BlockSize*1024;
                    % fseek(fid, fid_position);
                else           
                   
                    disp([kk,jj,llDiff,ii]);
                    tempBufferRare(1+posDataRead:readSize*2+posDataRead) = fread(fid, readSize*2, 'int32', dataFormat);  
                    posDataRead = posDataRead + readSize*2;
                   % position_fid = ftell(fid);
                   % display(sprintf('the current file position is %d  \n', position_fid));
                    
                    if flagKFormat == true; 
                       garbage = fread(fid, skipSize, 'int32', dataFormat); 
                    end;
%                      fid_position = ftell(fid);
%                     display(sprintf('the current file position after skip garbage is %d  \n', fid_position));
                end;   
                 
            end;   % end of: for ii = 1:(rareFactor+navFactor)
            
            if (flagSusumu == 1)  &&  (navFactor>0)
                if (kk == 1)  && (jj==1)
                   phaseBaseline(1:length(phaseOdd), llDiff) = calcPhaseBaseline(phaseOdd , phaseEven); 
                 
                       
                end;

                phaseShift  =  calcPhaseShift(phaseOdd, phaseEven, phaseBaseline(1:length(phaseOdd), llDiff));   % phaseSHift is calculated for this segment
            end;
           % phaseShift(1:length(phaseOdd),
           % (kk-1)*segmentSize+floor(jj/segmentSize)+mod(jj,segmentSize), llDiff) = calcPhaseEvolution(phaseOdd , phaseEven); 
          
           if (navFactor>0) && (dummy_counter(llDiff) < 120)  
            phaseOdd1 = fft(phaseOdd);
            phaseEven1 = fft(phaseEven);
            phaseInfoOdd(dummy_counter(llDiff), llDiff) = angle( phaseOdd1(5));
            phaseInfoEven(dummy_counter(llDiff), llDiff) = angle( phaseEven1(5));
            dummy_counter(llDiff) = dummy_counter(llDiff) + 1;
           end;
        
           if (flagSusumu == 1)   && (navFactor>0)
                for iiCorrection = 1:rareFactor
                    if (mod(iiCorrection,2)== 1)
                        tempBufferRare( 1+(iiCorrection-1)*readSize*2:readSize*2*iiCorrection ) = ...
                          correctPhaseShift( tempBufferRare( 1+(iiCorrection-1)*readSize*2:readSize*2*iiCorrection ),  phaseShift(1:length(phaseOdd),1), jj, kk);
                    else
                        tempBufferRare( 1+(iiCorrection-1)*readSize*2:readSize*2*iiCorrection ) = ...
                          correctPhaseShift( tempBufferRare( 1+(iiCorrection-1)*readSize*2:readSize*2*iiCorrection ),  phaseShift(1:length(phaseOdd),2), jj, kk);
                    end;
                end;
           end;

            
            count = fwrite(fid_rare, tempBufferRare,  'int32', dataFormat);
            if count ~= readSize*2*rareFactor
                error(' writing error into rare fid file ');
                return;
            end;       
        end;  % end of llDiff = 1:diffSize
 end;   %  jj = 1:(rareFactor+navFactor):(phaseSize*(rareFactor+navFactor)/rareFactor)
end;

if (flagPlotPhase == 1)  && (navFactor>0)
    for llDiff = 1:size(phaseInfoOdd, 2)
        plot(1:length(phaseInfoOdd), phaseInfoOdd(:,llDiff), '*r', 'MarkerSize', 10);
        hold on;
        plot(1:length(phaseInfoEven), phaseInfoEven(:, llDiff), '*b', 'MarkerSize', 10);
        hold on;
        plot(1:length(phaseInfoEven), (phaseInfoEven(:,llDiff)+phaseInfoOdd(:,llDiff))/2, 'g', 'LineWidth', 2);
        title(' the phase evolution of the navigator echoes of the consecutive segments ');
        % text('odd echoes');
        % text('even echoes');
        title(sprintf(' phase envolution of the %d diffSize ', llDiff));
        hold off;   
        pause(10);
    end;
end;

clear llDiff   % kk, ii, jj;
 clear kk ii jj;
 
 fclose(fid);
 fclose(fid_nav);
 fclose(fid_rare);
 
 
 
 %--------------  end of main function ----------------------------
 
function [ph0  ph1] = calcPh0Ph1( bufferNav, first, ii)
% this function is to calculate the zero order/ first order of phase for
% navigator echo

numFit = 5;
numFit = floor(numFit/2);

firstNav = first(1:2:end) + j*first(2:2:end);
signalMax1 = max(abs(firstNav(:)));
ttFirstNav = find(abs(firstNav(:)) == signalMax1);

signalNav = bufferNav(1:2:end) + j*bufferNav(2:2:end);
signalMax2 = max(abs(signalNav(:)));
ttSignalNav = find(abs(signalNav(:))==signalMax2);   % || find(abs(signalNav(:)<signalMax2/10)));

% ttDist =  ttSignalNav - ttFirstNav;
% firstNav = circshift(firstNav, [0 ttDist]); 

phaseDiff = angle(signalNav(ttSignalNav-numFit:ttSignalNav+numFit) - angle(firstNav(ttSignalNav-numFit:ttSignalNav+numFit)));
%phaseDiff = phase(phaseDiff);
% anglePlot(signalNav, ii);

% arrayNav = (1:length(signalNav)).'; 
% phaseDiff(ttSignalNav) = [];
arrayNav  = ttSignalNav-numFit:ttSignalNav+numFit;  
arrayNav = [ones(length(arrayNav), 1)  (arrayNav).' ];
phaseDiff = (phaseDiff).';
phaseDiff = phaseUnwrap(phaseDiff);
[ph0  ph1] = weighLinFit(phaseDiff, signalNav(ttSignalNav-numFit:ttSignalNav+numFit),ttSignalNav, numFit );
% b = regress(phaseDiff, arrayNav, 0.1);
% ph0 = b(1);
% ph1 = b(2);    % unit is cm
%*********    end of function calcPh0Ph1  *************

function dataCorrect = correctPh0Ph1(dataInput,  ph0,  ph1,jj,kk)
dataInputTemp = dataInput(1:2:end) + j*dataInput(2:2:end);
arrayData = 1:length(dataInputTemp);   %linspace(-0.5, 0.5, length(dataInputTemp)); 
dataInputTemp = dataInputTemp.*exp(-1*j*(ph0+arrayData*ph1));  % +arrayData*ph1
dataCorrect(1:length(dataInput)) = 0;    %zeros(length(dataInput),1);
dataCorrect(1:2:end) = real(dataInputTemp(:));
dataCorrect(2:2:end) = imag(dataInputTemp(:));

%*********    end of function correctPh0Ph1  *************

function anglePlot(signalNav, ii)
signalMax = max(abs(signalNav(:)));
tt = find(abs(signalNav(:))< signalMax/20);
signalNav(tt) = 0;
plot(abs(signalNav));
title(sprintf('the %d navigator echo ',  mod(ii,2)));
pause(1.5);
%**************** end of function of anglePlot **********************


function outputPhase = phaseUnwrap(inputPhase)

count1 = 1;
count2 = 1;
count3 = 1;
count4 = 1;
outputPhase = inputPhase;
%outputPhase(ceil(length(inputPhase)/2)) = inputPhase(ceil(length(inputPhase)/2));

for ii = ceil(length(inputPhase)/2)+1:length(inputPhase)
   % outputPhase(ii) = inputPhase(ii);
    if (inputPhase(ii) - inputPhase(ii-1))>pi
        for jj1 = ii:length(inputPhase)
            outputPhase(jj1) = inputPhase(jj1)-2*pi*count1;
        end;
         count1 = count1 + 1;
    end;
    
    if (inputPhase(ii) - inputPhase(ii-1))< -1*pi
        for jj1 = ii:length(inputPhase)
            outputPhase(jj1) = inputPhase(jj1)+ 2*pi*count2;
        end;
        count2 = count2+1; 
    end;
end;

for ii = ceil(length(inputPhase)/2)-1:-1:1
    % outputPhase(ii) = inputPhase(ii);
    if (inputPhase(ii) - inputPhase(ii+1))>pi
        for jj1 = ii:-1:1
            outputPhase(jj1) = inputPhase(jj1)-2*pi*count3;
        end;
         count3 = count3 + 1;
    end;
    
    if (inputPhase(ii) - inputPhase(ii+1))< -1*pi
        for jj1 = ii:-1:1
            outputPhase(jj1) = inputPhase(jj1)+ 2*pi*count4;
        end;
        count4 = count4+1; 
    end;
end;
%******************  end of function of phaseUnwrap  ********************
   
function [ph0  ph1] = weighLinFit(phaseDiff, weighed, ttSignalNav, numFit )
if length(phaseDiff) ~= length(weighed)
    error(' there is error in the length of input data');
    return;
end;

wFactor = sqrt(abs(weighed));
wFactor = floor(wFactor/min(wFactor(:)));

numIndex = (ttSignalNav-numFit):(ttSignalNav+numFit);


            
arrayData = [];
posIndex = 0;
for ii = 1:length(phaseDiff)
    temp = repmat(phaseDiff(ii), [wFactor(ii) 1]);
    arrayData = [arrayData;temp];
    arrayIndex(1+posIndex:posIndex+wFactor(ii)) = numIndex(ii);
    posIndex = posIndex + wFactor(ii);
end;
arrayIndex = arrayIndex.';
arrayIndex = [ones(length(arrayIndex),1)  arrayIndex];
x = arrayIndex\arrayData;


% numIndex = [ones(length(numIndex),1)   numIndex.'];
% x1 = numIndex\phaseDiff;

ph0 = x(1);
ph1 = x(2);


% plot(numIndex, phaseDiff, '*r', ...
%                 'LineWidth',0.00001,...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','g', ...
%                  'MarkerSize',10);
% hold on;
% arrayFit = ttSignalNav-numFit-3:0.5:ttSignalNav+numFit+3;
% plot(arrayFit, ph0+arrayFit*ph1, 'r');
% axis([ttSignalNav-numFit-4  ttSignalNav+numFit+4  -2*pi  2*pi]);
% title(' curve fitting for phase evolution ');
% hold off;
% 
% ph0;
              
%***************  end of function  weighLinFit  *****************

function  phaseShift = calcOddEvenPhase(phaseOdd , phaseEven)

if length(phaseOdd) ~= length(phaseEven)
    error(' the length is not the same');
    return;
end;
profileOdd = fft(phaseOdd);
profileEven = fft(phaseEven);
%tt = find(abs(profileOdd(:)) == max(abs(profileOdd(:))));
phaseShift = angle(profileOdd(:)) - angle(profileEven(:));    % uniy is the radian
phaseShift = phaseShift/2;
%**************  end of the function of calcOddEvenPhase ***************

function  phaseShift = calcPhaseShift(phaseOdd , phaseEven, phaseBaseline)

if length(phaseOdd) ~= length(phaseEven)
    error(' the length is not the same');
    return;
end;
profileOdd = fft(phaseOdd);
profileEven = fft(phaseEven);
%tt = find(abs(profileOdd(:)) == max(abs(profileOdd(:))));
phaseShift(1:length(phaseOdd),1) = phaseBaseline - angle(profileOdd(:));    % uniy is the radian
phaseShift(1:length(phaseOdd),2) = phaseBaseline - angle(profileEven(:));    % uniy is the radian
% tt = find(abs(profileOdd(:))< 0.2*max(abs(profileOdd(:))));
% maxProfileOdd = max(abs(profileOdd(:)));
% maxProfileEven = max(abs(profileEven(:)));
% profileOdd = profileOdd*(maxProfileEven/maxProfileOdd);
% phaseShift(1:length(phaseOdd),1) = (angle(profileOdd.*conj(profileEven)))/2;   %(angle(profileOdd) - angle(profileEven))/2;   % 
% phaseShift(1:length(phaseOdd),2) = -1*(angle(profileOdd.*conj(profileEven)))/2;   %-1*(angle(profileOdd) - angle(profileEven))/2;  %(angle(profileOdd.*conj(profileEven)))/2;    % 
% phaseShift(tt,1) = 0;
% phaseShift(tt,2) = 0;

%**************  end of the function of calcPhaseEvolution ***************

function outputSignal =  correctPhaseShift( inputSignal,  phaseShift, jj, kk)
outputSignal(1:length(inputSignal)) = 0; 
tempSignal = inputSignal(1:2:end) + j*inputSignal(2:2:end);
tempSignal = fft(tempSignal);
phaseShift = phaseShift.';
tempSignal = tempSignal.*exp(j*phaseShift); 
% tempSignal = ifft(tempSignal);
outputSignal(1:2:end) = real(tempSignal(:));
outputSignal(2:2:end) = imag(tempSignal(:));
%**************  end of the function of calcOddEvenPhase ***************

function    phaseBaseline = calcPhaseBaseline(phaseOdd , phaseEven)

if length(phaseOdd) ~= length(phaseEven)
    error(' the length is not the same');
    return;
end;
profileOdd = fft(phaseOdd);
profileEven = fft(phaseEven);
%tt = find(abs(profileOdd(:)) == max(abs(profileOdd(:))));
phaseBaseline = angle(profileOdd(:)) + angle(profileEven(:));    % uniy is the radian
phaseBaseline = phaseBaseline/2;

