function acquisition_time_dicom(flist)

cur_d = pwd;

for i=1:length(flist)
    
    cd(fullfile(cur_d,flist{i}));
    dir_s = dir('*.ima');
    cmd = sprintf('dicom_hdr %s | grep -e "Acquisition Time"',dir_s(1).name);
    [e,s]=unix(cmd);
    disp(s);
end
cd(cur_d);