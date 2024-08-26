function average_pcscans(fpat,scanIDs)
% mask_PAinPC(fpat,scanIDs,wm_mask,interp_factor)
warning off;
nscans=length(scanIDs);

for i=1:nscans
    dname=sprintf(fpat,scanIDs(i)); 
    data(:,:,:,:,i)=ri(dname,1);    
end
tmp=load(dname);

tmp.d=mean(data,5);

 fpat2=strrep(fpat,'%d','%s');    
 fpat2=strrep(fpat2,'%02d','%s');
 fpat2=strrep(fpat2,'%03d','%s');
 fpat2=strrep(fpat2,'%04d','%s');
    
savename=sprintf(fpat2,[num2str(scanIDs(1)),'_',num2str(scanIDs(end))]);

    
save(savename,'-struct','tmp');



    
    
    
