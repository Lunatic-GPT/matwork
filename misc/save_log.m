function save_log(varargin)

fid=fopen(sprintf('log_%s.txt',varargin{1}),'w');
    fprintf(fid,'%s(',varargin{1}); 
for i=2:length(varargin)
    
    fprintf(fid,'%s,',num2str(varargin{i})); 
end

    fprintf(fid,')'); 
    fclose(fid);
