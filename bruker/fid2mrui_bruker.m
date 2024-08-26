function fid2mrui_bruker(dname,prefix)


if iscell(dname)
    
    
    
    for i=1:length(dname)
       if i==1
        d=readb_fid(dname{i});
        
        sw=readbPar(fullfile(dname{i},'method'),'PVM_SpecSWH',true);

       else
           
           tmp=readb_fid(fullfile(dname{i}));
           d=cat(2,squeeze(d),squeeze(tmp));
       end
       
      
        
    end
else 
    
    
d=readb_fid(dname);

sw=readbPar(fullfile(dname,'method'),'PVM_SpecSWH',true);


end


write_mrui2(d,prefix,sw);
