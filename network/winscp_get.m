function status=winscp_get(fname,server,dname,ntrials)
%  winscp_get(killdevil,dname,fname)
if ~exist('ntrials','var')
    ntrials = Inf;
end

if ~exist('server','var')
    server='longleaf';
end

if ~exist('dname','var')
   dname= '/nas/longleaf/home/zongx/pine/temp';
end

cmd=sprintf('%s   "cd %s"',cmd_connect(server),dname);

cmd=sprintf('%s "lcd ""%s""" "get %s"',cmd,strrep(pwd,'\','/'),fname);
[tmp,status]=fail_free(cmd,ntrials);


    
    
    
    