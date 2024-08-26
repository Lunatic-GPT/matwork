function write_DiffusionVectors(vec,fname)
% write_DiffusionVectors(vec,suffix)
prefix=strtok(fname,'.');
fid=fopen(sprintf('%s.dvs',prefix),'w');

if iscell(vec)
    
    for j=1:length(vec)
        fprintf(fid,'[directions=%d]\n',size(vec{j},1));
        fprintf(fid,'coordinatesystem=xyz\n');
        fprintf(fid,'normalisation = none\n');
        
        
        for i=1:size(vec{j},1)
            
            fprintf(fid,'vector[%d]  = (%f, %f, %f)\n',i-1,vec{j}(i,1),vec{j}(i,2),vec{j}(i,3));
        end
        
        if j<length(vec)
            fprintf(fid,'\n');
            fprintf(fid,'\n');
            fprintf(fid,'\n');
            fprintf(fid,'\n');
        end
    end
    
    
else
    
    fprintf(fid,'[directions=%d]\n',size(vec,1));
    fprintf(fid,'coordinatesystem=xyz\n');
    fprintf(fid,'normalisation = none\n');
    
    
    for i=1:size(vec,1)
        
        fprintf(fid,'Vector[%d]  = (%f, %f, %f)\n',i-1,vec(i,1),vec(i,2),vec(i,3));
    end
end

fclose(fid);

read_diffvectors(sprintf('%s.dvs',prefix)); % check output file;