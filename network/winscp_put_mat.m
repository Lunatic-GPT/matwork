function winscp_put_mat(m_name)
% put 
% fname: matlab file 

fname=which(m_name);
[dname,fname,ext]=fileparts(fname);
fname=[fname,ext];

ind=strfind(dname,'matwork');
cmd='winscp.exe /command "open longleaf" "option batch on"';
dname=strrep(dname,'\','/');
%cmd=sprintf('%s   "cd /nas2/home/z/o/zongx/matwork/%s"',cmd,dname(ind+8:end));
cmd=sprintf('%s   "cd /nas/longleaf/home/zongx/matwork/%s"',cmd,dname(ind+8:end));

cmd=sprintf('%s "lcd ""%s""" "put %s"',cmd,dname,fname);
fail_free(cmd);


function results=fail_free(cmd)
    status=1;
    while status~=0
       fprintf('.');
       [status,results]=system(cmd);
       fprintf('+\n');
    end
    
    
    
    