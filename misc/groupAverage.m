function [a,a_std] = groupAverage(flist,var)
 
nsub = length(flist);
 

 for i=1:nsub
    load(flist{i},var);
    
    if i==1
      arr = eval(var);
      nd = length(size(arr));
    else
        arr = cat(nd+1,arr,eval(var));
    end
    
 end
 
 a = mean(arr,nd+1);
 a_std = std(arr,0,nd+1)/sqrt(nsub);
 
 
 