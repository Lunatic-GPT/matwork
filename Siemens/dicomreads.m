function [c,info] = dicomreads(fname)
% c = dicomreads(fname); read spectroscopy data 
 
if ~isdir(fname)
fd=dicom_open(fname);

c=dicom_get_spectrum_siemens(fd);
fclose(fd);

info=dicomreads_header(fname);

else
    
    str=dir([fname,'\*.IMA']);
    if isempty(str)
     str=dir([fname,'\*.dcm']);
    end    
    
    for i=1:length(str)
    fd=dicom_open([fname,'\',str(i).name]);

    c(:,i)=dicom_get_spectrum_siemens(fd);
    fclose(fd);

    end
info=dicomreads_header([fname,'\',str(i).name]);


end

