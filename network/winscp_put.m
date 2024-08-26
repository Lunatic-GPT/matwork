function winscp_put(src_name,server,dest_dir)
% winscp_put(src_name,server,dest_dir) 


if ~exist('server','var')
   server='longleaf'; 
end

if ~exist('dest_dir','var')
   dest_dir= '/nas/longleaf/home/zongx/pine/temp';
end

if isempty(fileparts(src_name))
 src_name=fullfile(pwd,src_name);
end

dname=fileparts(src_name);

dname=strrep(dname,'\','/');
fname=filename(src_name);

cmd=sprintf('%s "cd %s" "lcd ""%s""" "put ""%s""" "exit"',cmd_connect(server),dest_dir,dname,fname);
fail_free(cmd);




    
    
    