function parfile_convert(pfile,parray,fout)

parray = str2cell(parray);

fid = fopen(fout,'a+');

a=textscan(fid,'%s');

for i=1:length(parray)
    ind = strmatch(parray{i},a{1},'exact');
    if ~isempty(ind)
        fprintf('%s exist, skipped\n',parray{i});
        continue;
    end
    
    v = readPar(fileparts(pfile),parray{i});
    
    
    if isa(v,'char')
        
     fprintf(fid,'%s %s\n',parray{i},v(2:end-1));
    
    else
     fprintf(fid,'%s %f\n',parray{i},v);
    end
end



