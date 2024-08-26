function dicom_keepone(pat)


a=find_pattern(pat);

for i=1:length(a)
    
    if ~exist(a{i},'dir')
         continue;
    end
    
   do_keepone_folder(a{i}); 
    

  
    
end



function do_keepone_folder(dname)

 b=dir2(dname,2);
      
 for i=2:length(b)    
    delete(fullfile(b(i).folder,b(i).name));
    fprintf('delete %s %s\n',b(i).folder,b(i).name);
 end
  
 c=dir2(dname,1);
 for j=1:length(c)
      do_keepone_folder(fullfile(c(j).folder,c(j).name)); 
 end
   
 
 
 