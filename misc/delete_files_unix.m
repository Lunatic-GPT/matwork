function delete_files_unix(pattern,d_name)
%delete_files(pattern,d_name)
if ~exist('d_name','var');
    d_name = pwd;
end

cmd = sprintf('find %s -name "%s"',d_name,pattern);
[status,results]=unix(cmd);
str = textscan(results,'%s','delimiter','\n');

for i=1:length(str{1})
    cmd = sprintf('rm -rf ''%s''',str{1}{i});
    unix(cmd);
    disp(str{1}{i});
    disp(i);
end