function dod(pattern)

if ~exist('pattern','var')
    pattern='';
end

cmd=sprintf('dir /od %s',pattern);

system(cmd);