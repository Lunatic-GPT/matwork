cmd = 'find . -name .svn';
[status,results]=unix(cmd);
str = textscan(results,'%s','delimiter','\n');

for i=1:length(str{1})
    cmd = sprintf('rm -rf ''%s''',str{1}{i});
    unix(cmd);
    disp(str{1}{i});
    disp(i);
end


cmd = 'find . -name "*.bmp"';
[status,results]=unix(cmd);
str = textscan(results,'%s','delimiter','\n');

for i=1:length(str{1})
    cmd = sprintf('mv ''%s'' ../Migraine/Stimulus_Images',str{1}{i});
    unix(cmd);
    disp(str{1}{i});
    disp(i);
end