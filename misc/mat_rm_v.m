function mat_rm_v(mat,vname)


tmp=load(mat);

try
tmp=rmfield(tmp,vname);
    

save(mat,'-struct','tmp');

catch
   warning('Field %s does not exist',vname); 
end