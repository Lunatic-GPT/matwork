function [a,a_std] = groupAverage(dlist,fname,var)
 
nsub = length(dlist);
 
 for i=1:nsub
    load(fullfile(dlist{i},fname),var);
    
    if i==1
        arr = zeros([size(eval(var)),nsub]);
    end
    
    arr(:,:,i) = eval(var);
    
 end
 
 a = mean(arr,3);
 a_std = std(arr,0,3);
 
 
 