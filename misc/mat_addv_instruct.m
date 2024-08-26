function mat_addv_instruct(mat,str)


tmp=load(mat);


name=fieldnames(str);

for i=1:length(name)
    val=getfield(str,name{i});
 eval(sprintf('tmp.%s=val;',name{i}));
end

save(mat,'-struct','tmp');