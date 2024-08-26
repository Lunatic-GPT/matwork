    function cmd=cmd_connect(server)
    % dname: file or dir name; relative to netscr (killdevil) and pine (longleaf)

         cmd=sprintf('winscp.exe /command "open %s" "option batch on"',server); 
    
        