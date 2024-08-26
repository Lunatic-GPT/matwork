function mat2dcm(matname,dcmname)
%mat2dcm(matname,dcmname)
tic;


dir_str=dir([dcmname,'/*.dcm']);
if length(dir_str)==0
dir_str=dir([dcmname,'/*.IMA']);
end    

dinfo=dicominfo(fullfile(dcmname,dir_str(1).name));

b=ri(matname);
b=(b-min(b(:)))*4095/(max(b(:))-min(b(:)));

prefix=strtok(matname,'.');
mkdir(prefix);
dicomwrite(uint16(b),sprintf('%s.dcm',prefix),dinfo);
movefile( sprintf('%s.dcm',prefix),prefix);
 
%% save the average time series for each trial type.    

    disp([mfilename ' finish in ', num2str(toc), ' s']);