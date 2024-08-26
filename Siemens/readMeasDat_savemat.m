function [Data,navData]=readMeasDat_savemat(fname,nline)
% VB: [Data,PMUTimeStamp,freePara,navData]=readMeasDat_savemat(fname,nline)
% VE:  [Data,navData]=readMeasDat_savemat(fname)
% Same as readMeasDat but save the extracted data as a mat file in the same directory before
% returning the data.  
% If the mat file already exists, data will be retrieved from the mat file
% instead for speed.
% If % Navigator (mdhBitFields.MDH_REFPHASESTABSCAN || mdhBitFields.MDH_RTFEEDBACK)
% data are available, they will be saved separately.

if ~exist('nline','var')
    nline=Inf;
end

prefix=strtok2(fname,'.');
if ~exist([prefix,'.mat'],'file')
    
    if strcmp(rawDataVersion(fname),'vb')
        [Data,navData] = readMeasDat(fname,nline,0,true);        
       
    else
        [Data,navData] = readMeasDatVE11(fname,true);
        
    end
    
    save([prefix,'.mat'],'-struct','Data','-v7.3');  
    if ~isempty(navData)        
      save([prefix,'_FatNav.mat'],'-struct','navData','-v7.3');
    end
 
else
    if nargout>=1
     Data=load([prefix,'.mat']);
    end
    
    if exist([prefix,'_FatNav.mat'],'file') && nargout>1
        navData=load([prefix,'_FatNav.mat']);
    end 
end
