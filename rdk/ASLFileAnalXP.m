function ASLFileAnalXP(eyedatafile, nTR)
%function ASLAnalXP()

echo on
%Xinmiao, modifying from c:\Program Files\ASL Eye Tracker 6000\SDK\Matlab\ASLFileSample.m,
%05/11/09

%Params for calibration, specific for scanner #1 setup
hfactor = (atan(23/78)/pi*180)/(215-42);
vfactor = (atan(18/78)/pi*180)/(193-46);
centerX = 125;
centerY = 122;

TRTriggerXdat = 132;

svr = actxserver('ASLFileLib2.ASLFile')

% uncomment the next line to see all methods available from the DLL
%svr.methodsview

% Open file
[segmentCount, itemNames] = invoke(svr, 'OpenFileForReading', eyedatafile)
%[segmentCount, itemNames] = invoke(svr, 'OpenFileForReading', 'AIONMSRFL.eyd')

tAll = {};
hposAll = {};
vposAll = {};
xdatAll = {};
for iSeg = 1:segmentCount
    segSize = svr.GetSegmentInfo(iSeg);
    record = svr.ReadDataRecord;
    tsec0 = convertToSec(record{3});
    t = 0;
    hpos = 0;
    vpos = 0;
    xdat = 0;
    for iRec = 2:segSize
        record = svr.ReadDataRecord;
        t(iRec-1) = convertToSec(record{3}) - tsec0;
        hpos(iRec-1) = record{12};
        vpos(iRec-1) = record{13};
        xdat(iRec-1) = record{10};
    end

    indLastTR = t(xdat == TRTriggerXdat);
    indLastTR = indLastTR(length(indLastTR)-5); % 6 fields for one TR trigger
    hpos = (hpos(t<=indLastTR) - centerX)*hfactor;
    vpos = (vpos(t<=indLastTR) - centerY)*vfactor;
    xdat = xdat(t<=indLastTR);
    t = t(t<=indLastTR);
    t = t - indLastTR + nTR*2 -2;

    tAll{iSeg} = t;
    hposAll{iSeg} = hpos;
    vposAll{iSeg} = vpos;
    xdatAll{iSeg} = xdat;
end
    

% Close file
svr.CloseFile

% Release the memory
svr.delete

eval(['save ' eyedatafile(1:length(eyedatafile)-4) ' *All']); 


function tsec = convertToSec(rec3);
% Converting the 3rd item of record time to sec.

h = str2num(rec3(1:2));
m = str2num(rec3(4:5));
s = str2num(rec3(7:12));
tsec = h*3600+m*60+s;
