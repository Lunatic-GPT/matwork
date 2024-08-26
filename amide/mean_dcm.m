function mean_dcm(dir_name)
% generate the mean of all dicom files in the dir_name

a=ri(dir_name,1);

a2=mean(a,4);

dir_str=dir(dir_name);

in=dicominfo(fullfile(dir_name,dir_str(end).name));

mkdir([dir_name,'_mean']);

dicomwrite(uint16(a2),[dir_name,'_mean.IMA'],in);

movefile([dir_name,'_mean.IMA'],[dir_name,'_mean\']);