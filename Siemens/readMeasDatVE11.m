function [sdata,snavData] = readMeasDatVE11(filename,remove_os)
%[actualData,ushulPMUTimeStamp,ushPartition,ushSlice] = readMeasDat(filename,max_nlines,offset[,remove_os])
% filename: data file name
% sdata and snavData has the following fields:
% Line, Partition, Set, Slice, Channel,  Repetition, Data;

% Line and Partition are 0 based;


if ~exist('remove_os','var')
    remove_os=false;
end
d=mapVBVD(filename);

sdata=extractTwixObj(d{2}.image);

if isfield(d{2},'refscan')
   srefData=extractTwixObj(d{2}.refscan);    
   
   lin_par_ref=[srefData.Line(:),srefData.Partition(:)];
   lin_par=[sdata.Line(:),sdata.Partition(:)];
   
   [~,id]=setdiff(lin_par_ref,lin_par,'rows');
   sdata.Line=cat(2,sdata.Line,srefData.Line(id));
   sdata.Partition=cat(2,sdata.Partition,srefData.Partition(id));
   sdata.Slice=cat(2,sdata.Slice,srefData.Slice(id));
   sdata.Data=cat(3,sdata.Data,srefData.Data(:,:,id));
   sdata.Repetition=cat(2,sdata.Repetition,srefData.Repetition(id));
   sdata.Set=cat(2,sdata.Set,srefData.Set(id));
   sdata.iceParam=cat(2,sdata.iceParam,srefData.iceParam(:,id));
   sdata.freeParam=cat(2,sdata.freeParam,srefData.freeParam(:,id));
   sdata.timestamp=cat(2,sdata.timestamp,srefData.timestamp(id));
   
end

if isfield(d{2},'RTfeedback')
   snavData=extractTwixObj(d{2}.RTfeedback);    
else
   snavData=[]; 
end




   if remove_os
            tmp=ifft1c(sdata.Data,1);
            len=size(tmp,1);
            sdata.Data=tmp(len/4+1:end-len/4,:,:);    
            
            if ~isempty(snavData)
              tmp=ifft1c(snavData.Data,1);
              len=size(tmp,1);
              snavData.Data=tmp(len/4+1:end-len/4,:,:);    
            end
            
             if ~isempty(srefData)
              tmp=ifft1c(srefData.Data,1);
              len=size(tmp,1);
              srefData.Data=tmp(len/4+1:end-len/4,:,:);    
            end
   end




