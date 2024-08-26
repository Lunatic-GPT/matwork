function run_cmd_remote(cmd,server,dname,ntrials)
% dname: relative folder name in /netscr or /pine; will be created if not
% exist on server
% fname: name of the shell script in the local directory.



if ~exist('server','var') || isempty(server)
    server='longleaf';
end

if ~exist('dname','var') || isempty(dname)
   dname = '/nas/longleaf/home/zongx/pine/temp';
end

cmd=str2cell(cmd);
cmd_str='';
for i=1:length(cmd)
    cmd_str=sprintf('%s "call %s"',cmd_str,cmd{i});
end

cmd_all=sprintf('%s "cd %s" %s "exit"',cmd_connect(server),dname,cmd_str);

fail_free(cmd_all,ntrials);


    
    
    
        