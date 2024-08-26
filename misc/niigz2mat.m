function niigz2mat

dir_str=dir('*.gz');

for i=1:length(dir_str)
    
load_untouch_niigz(dir_str(i).name);

end