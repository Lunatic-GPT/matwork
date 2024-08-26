function mat_addv(mat,vname,val)


tmp=load(mat);

eval(sprintf('tmp.%s=val;',vname));
    

save(mat,'-struct','tmp');