function bsub_FatNav(fname,lstart,lend,lstep)

mid=strtok_no(fname,'_',2);
fid=fopen(unique_name('bscript.sh'),'w');
if exist('lstart','var')
    
    
    ind=lstart:lstep:lend;
    ind2=ind+lstep-1;
    if ind2(end)>lend
        ind2(end)=lend;
    end
    
    for i=1:length(ind)
        fprintf(fid,'bsub -M 16 -o FatNavReconLog.%%J matbgk "do_recon_FatNav(''%s'',%d:%d)" FatNavRecon%s_%d.log\n',fname,ind(i),ind2(i),mid,i);
        
    end
    
    fprintf(fid,'bsub -M 16 matbgk "do_volreg_afterFatNav(''%s'',%d:%d:%d)" FatNavRecon%s.log',mid,lstart,lstep,lend,mid);
    
    
else
    fprintf(fid,'bsub -M 16 -o FatNavReconLog.%%J matbgk "do_recon_FatNav(''%s'')" FatNavRecon%s.log\n',fname,mid);
    
    fprintf(fid,'bsub -M 16 matbgk "do_volreg_afterFatNav(''%s'',1)" FatNavRecon%s.log',mid,mid);
    
end

fclose(fid);