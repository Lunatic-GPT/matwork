function ASLFileAnal(eyedatafile)
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

svr = actxserver('ASLFileLib2.ASLFile');

[segmentCount, itemNames] = invoke(svr, 'OpenFileForReading', eyedatafile);

record_arr = cell(1,segmentCount);
time_arr = cell(1,segmentCount);
for iSeg = 1:segmentCount
    disp(iSeg);
    segSize = svr.GetSegmentInfo(iSeg);
    record = zeros(20,segSize);
    time = zeros(12,segSize,'uint8');
    for iRec = 1:segSize
        tmp = svr.ReadDataRecord;
        iit_arr = [1,2,4:21];
        for iit = 1:20
          record(iit,iRec) = tmp{iit_arr(iit)};
        end
        time(:,iRec) = tmp{3};
    end
record_arr{iSeg} = record;
time_arr{iSeg} = time;
end
    
svr.CloseFile

% Release the memory
svr.delete

prefix = strtok(eyedatafile,'.');
save(prefix,'record_arr','itemNames','time_arr'); 

