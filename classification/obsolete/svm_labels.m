function svm_labels(seq,fname_label)
fid = fopen(fname_label,'w');


n = length(seq);

for i=1:n
   if seq(i) == 1 
    fprintf(fid,'1\n');
   elseif seq(i) == 2
    fprintf(fid,'-1\n');
   else 
     fprintf(fid,'9999\n');  
   end
end

fclose(fid);
