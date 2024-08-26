function mkdir_remote(dname,server)
% dname: relative folder name in /netscr or /pine; will be created if not
% exist on server
% fname: name of the shell script in the local directory.

cmd=cmd_connect(server);

if ~dir_exist(dname,server)
    cmd2=sprintf('%s "mkdir %s"',cmd,dname);
    fail_free(cmd2);
end


    
    
    
        