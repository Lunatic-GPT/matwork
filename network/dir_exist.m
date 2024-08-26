    function res=dir_exist(dname,server)
    % dname: file or dir name; relative to pine (longleaf)
    [dname,fname,ext]=fileparts(dname);
    
    if isempty(dname)||dname(1)~='/'
        dname = ['/pine/scr/z/o/zongx/',dname];
    end
    
    cmd=sprintf('winscp.exe /command "open %s" "option batch on" "cd %s"',server,dname); 
    
    fname=[fname,ext];
    
        cmd2=[cmd,' "call ls>winscp_temp.txt"'];
        
        results=fail_free(cmd2);
        
          
        cmd2=sprintf('%s "lcd ""%s""" "get winscp_temp.txt"',cmd,strrep(pwd,'\','/'));
        fail_free(cmd2);
        
        fid=fopen( 'winscp_temp.txt','r');
        a=textscan(fid,'%s');
        fclose(fid);
        
        if ~isempty(strmatch(fname,a{1},'exact'))        
            res=true;
        else
            res=false;
        end
        