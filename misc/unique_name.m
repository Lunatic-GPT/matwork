function res=unique_name(fname)

[prefix,suf]=strtok2(fname,'.');

count='a'-1;
prefix2=prefix;

while 1
    if exist([prefix2,suf],'file')
        count=count+1;
        prefix2 =[prefix,'_',count];
    else
        break;
    end
    
    
end
res=[prefix2,suf];

