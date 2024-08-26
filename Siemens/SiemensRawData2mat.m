function SiemensRawData2mat(fpattern,nline)
% a wrapper for readMeasDat_savemat so that multiple files can be processed
% in one functional call.
% old name: extractFatNavData;

if ~exist('nline','var')
    nline=Inf;
end

fname_list=name4pat(fpattern);

for i=1:length(fname_list)
    fname=fname_list{i};   
    prefix=strtok2(fname,'.');
    if ~exist([prefix,'.mat'],'file')
        readMeasDat_savemat(fname,nline);       
    end
end


