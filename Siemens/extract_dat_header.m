function extract_dat_header(fname)


fid=fopen(fname,'rb');

prefix=strtok(fname,'.');

fid2=fopen([prefix,'_header.txt'],'w');
i=0;
while ~feof(fid)
    s=fgetl(fid);
    
    if length(s)>10 && ~contain_text(s)
        break;
    end
    
    fprintf(fid2,'%s\n',s);
   
end
fclose(fid2);
fclose(fid);




function res=contain_text(s)


if any(int32(s)<32 | int32(s)>126)
    
    if isempty(strfind(s,'XProtocol')) && isempty(strfind(s,'ASCCONV'))
       res=false;
    else
        res=true;
    end
    
else
    res=true;
end

    

