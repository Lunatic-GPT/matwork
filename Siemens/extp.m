function extp(fname)
%% extract protocol

if exist(fname,'dir')
    dir_str=dir(fname);
    
    i=0;
    
    prefix=fname;
    while 1
        if exist(fullfile(fname,dir_str(end-i).name),'dir')
            
            i=i+1;
            
        else
            fid=fopen(fullfile(fname,dir_str(end-i).name),'rb');
            break;
        end
    end
else
    fid=fopen(fname,'rb');
    prefix=strtok2(fname,'.');
end
%for i=1:length(s)
if exist([prefix,'.pro'],'file')
  return;
end
%   if strmatch(s{i},'ASCCONV')
fid2=fopen([prefix,'.pro'],'w');
start=0;
while ~feof(fid)
    s=fgetl(fid);
    if ~isempty(strfind(s,'ASCCONV'))
        start=start+1;
        
        if start==1
            fprintf(fid2,'### ASCCONV BEGIN ###\n');
            continue;
        else
            fprintf(fid2,'### ASCCONV END ###\n');
            break;
        end
        
    end
    
    if start==1
        
        fprintf(fid2,'%s\n',s);
        
    end
end

fclose(fid);
fclose(fid2);


