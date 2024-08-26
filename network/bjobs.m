function bjobs(killdevil)
%  bjobs(killdevil)

if isa(killdevil,'char')
    killdevil=str2num(killdevil);
end

if killdevil
run_cmd_remote('bjobs >bjobs.log','/netscr/zongx',1);

winscp_get(1,'/netscr/zongx','bjobs.log');
else
run_cmd_remote('sq >bjobs.log','/nas/longleaf/home/zongx/temp',0);

winscp_get(0,'/nas/longleaf/home/zongx/temp','bjobs.log');
    
end
!more bjobs.log

