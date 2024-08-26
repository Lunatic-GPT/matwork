function res=name4pat(pat,first_only,dir_only)
% name4pat(pat)
% if first_only is true, return a string, otherwise return a cell array
% even if it contains only one element

if ~any(pat=='*') %not a pattern
    res={pat};
    return;
end

if ~exist('dir_only','var')
    dir_only=false;
end

if ~exist('first_only','var')
    first_only=false;
end

dir_str=dir(pat);



if isempty(dir_str)
   res=[];
   return;
end


if length(dir_str)==1
 res=fullfile(dir_str(1).folder,dir_str(1).name);
 
 if dir_only
    if ~exist(res,'dir')
        res=[];
    end
 end
 
 if ~isempty(res)
     res={res};
 end
else
    res={};
    for i=1:length(dir_str)
        tmp=fullfile(dir_str(1).folder,dir_str(i).name);
        if dir_only
            if exist(tmp,'dir')
                res{end+1}=tmp;
            end
        else
            res{end+1}=tmp;
        end
    end
end

if first_only && ~isempty(res)
    res=res{1};
end
