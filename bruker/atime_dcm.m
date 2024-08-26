function atime_dcm(dname)

d=dir(fullfile(dname,'*.dcm'));

h=dicominfo(fullfile(dname,d(3).name));

fprintf('%s %s\n',h.AcquisitionDate,h.AcquisitionTime);

