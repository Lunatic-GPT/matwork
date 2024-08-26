function prob_prefix=pvsseg_remote(dpattern,server,action)
% action:
% 0: segment and download
% 1: download only
% 2: return prefix only
if ~exist('server','var') || isempty(server)
    server='andrew';
end

if ~exist('action','var')
    action=0;
end

dall=name4pat(dpattern);
dall=str2cell(dall);

for i=1:length(dall)
    fname=filename(dall{i});
    
    prefix=remove_suffix(fname,'.nii.gz');
    prob_prefix=['Prob_m2edn_nmsk_iter0_',prefix];
    if action==0
        if ~exist([prefix,'.nii.gz'],'file')
            dcm2nii(prefix);
        end
        dcmpar=filename_append(prefix,'DCMPar_');
        if exist(dcmpar,'file')
            a=load(dcmpar);
            if a.par.VoxelSize(1)<0.39  %interpolated and not yet undersampled
                undersample_image(sprintf('%s.nii.gz',prefix),[2,2,1]);
                movefile(sprintf('%s_ds.nii.gz',prefix),sprintf('%s.nii.gz',prefix));
            end
        end
        
        
        if exist([prob_prefix,'.nii.gz'],'file')
            continue;
        end
        disp(prefix);
        upload_data(prefix,server);    
        run_pvsseg(prefix,server);
         disp('Pause for 40 s'); 
       pause(40); % somehow the downloaded data was not completely, so wait 8 s.
    end
   
    if action<=1
     download_prob(prob_prefix,server);
    end
    
end

function upload_data(dname,server)

cmd=sprintf('winscp.exe /command "open %s" "option batch on" "cd /home/zongx/PVS_seg" "rm *.gz"',server);
cmd=sprintf('%s "lcd ""%s"""',cmd,pwd);

cmd=sprintf('%s "put ""%s.nii.gz"""',cmd,dname);
fail_free(cmd);




function run_pvsseg(dname,server)

cmd=sprintf('winscp.exe /command "open %s" "option batch on" "cd /home/zongx/PVS_seg"',server);
%cmd=sprintf('%s "lcd "%s""',cmd,pwd);
cmd=sprintf('%s "call pvsseg.py %s.nii.gz"',cmd,dname);
fail_free(cmd,1);


function download_prob(dname,server)

cmd=sprintf('winscp.exe /command "open %s" "option batch on" "cd /home/zongx/PVS_seg"',server);
cmd=sprintf('%s "lcd ""%s"""',cmd,pwd);
cmd=sprintf('%s "get %s.nii.gz"',cmd,dname);
fail_free(cmd);
