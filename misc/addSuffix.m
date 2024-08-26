function addSuffix(pat,suffix)

%rename_files(p,newp)
l=dir(['*',pat,'*']);

for i=1:length(l)
    
    name=l(i).name;
    
    movefile(name,[name,suffix]);
end
