function mat2text(matname,fieldname)

prefix=strtok2(matname,'.');
load(matname,fieldname);
y=eval(fieldname);

if ~exist(prefix,'dir')
    mkdir(prefix);
end
fid=fopen(fullfile(prefix,[fieldname,'.txt']),'w');

if isa(y,'cell')
    
    for i=1:length(y)
       
        for j=1:length(y{i})
           
           fprintf(fid,'%s ',num2str(y{i}(j)));  
            
        end
        
       fprintf(fid,'\n');
    end
    
else
   for i=1:size(y,1)
       for j=1:size(y,2)
           
           fprintf(fid,'%s ',num2str(y(i,j)));
           
       end
       fprintf(fid,'\n');
   end
    
end
fclose(fid);

