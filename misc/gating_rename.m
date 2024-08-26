function gating_rename(dur_min)
% rename gating output files with duration longer than dur_min seconds
% gating_rename(dur_min)
% dur_min: minimun duration of scans in s.  This will exclude files generated
% during prescans and reference scans.
dir_str = dir('RespData*');
dir_str2=dir('ECGData*');
iscan = 0;
for i=1:length(dir_str)
    a = load(dir_str(i).name);
    if length(a)>40*dur_min
        iscan = iscan+1;
        fprintf('Scan %d duration: %3.2f\n',iscan,length(a)/40);
        copyfile(dir_str(i).name,sprintf('Resp_s%d.1D',iscan));
        fprintf('Copy %s to %s\n',dir_str(i).name,sprintf('Resp_s%d.1D',iscan));
         copyfile(dir_str2(i).name,sprintf('ECG_s%d.1D',iscan));
        fprintf('Copy %s to %s\n',dir_str2(i).name,sprintf('ECG_s%d.1D',iscan));
    end
end
