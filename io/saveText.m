function saveText(mat,fname)
fid=fopen(fname,'w');

for i=1:size(mat,1)
    fprintf(fid,'%d   ',i);
    for j=1:size(mat,2)
   
      fprintf(fid,'%s  ',num2str(mat(i,j)));    
        
    end
    
    fprintf(fid,'\n');
end
fclose(fid);
