function res=TextLines_comments(fname,commentChar)

fid = fopen(fname);
tline = fgetl(fid);
count=0;
res={};
while ischar(tline)
    
    ind=find(tline==commentChar,1);
    if ~isempty(ind)
        tmp=strtrim(tline(1:ind-1));
    else
        tmp=strtrim(tline);
    end
    
    if ~isempty(tmp)
        count=count+1;      
        res{count}=tmp;
    end
    
    tline = fgetl(fid);
end
fclose(fid);