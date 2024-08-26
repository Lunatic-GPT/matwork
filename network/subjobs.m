function subjobs(fname,dname,killdevil)
% dname: relative folder name in /netscr or /pine; will be created if not
% exist on server
% fname: name of the shell script in the local directory.

cmd=cmd_connect(killdevil);

if ~dir_exist(dname,killdevil)
   mkdir_remote(dname,killdevil);
end

cmd=sprintf('%s "cd %s" "lcd %s" "put %s" "call source %s" "exit"',cmd,dname,pwd,fname,fname);
fail_free(cmd);


    
    
    
        