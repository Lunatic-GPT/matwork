function pvsseg_remote_longleaf(dpattern)

dall=name4pat(dpattern);
dall=str2cell(dall);

for i=1:length(dall)
    dname=dall{i};
   
       
    
    upload_data(dname);
    th=tic;
    run_pvsseg(filename(dname));

    
    download_prob(filename(dname));
    
end

function upload_data(dname)

cmd=sprintf('scp %s zongx@erwin "option batch on" "cd /home/zongx/PVS_seg"';
cmd=sprintf('%s "lcd ""%s"""',cmd,pwd);
cmd=sprintf('%s "put %s.nii.gz"',cmd,dname);
fail_free(cmd);




function run_pvsseg(dname)

cmd='winscp.exe /command "open erwin" "option batch on" "cd /home/zongx/PVS_seg"';
%cmd=sprintf('%s "lcd "%s""',cmd,pwd);
cmd=sprintf('%s "call pvsseg.py %s.nii.gz"',cmd,dname);
fail_free(cmd);


function download_prob(dname)

cmd='winscp.exe /command "open erwin" "option batch on" "cd /home/zongx/PVS_seg"';
cmd=sprintf('%s "lcd ""%s"""',cmd,pwd);
cmd=sprintf('%s "get prob_%s.nii.gz"',cmd,dname);
fail_free(cmd);
