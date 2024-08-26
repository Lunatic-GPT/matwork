function [i1,i2]=time_match(t1,t2,dtmax,max_t,min_t)
% [i1,i2]=time_match(t1,t2,dtmax,max_t,min_t)
% find the times in t1 and t2 that are within dtmax of each other and
% return their indices
% if a time matches with more than one time in the other array, then return
% the closest time in the other array.
% assume t1 and t2 are both sorted in the same order.
i1=[];
i2=[];
if length(t1)<length(t2)
 for i=1:length(t1)
    
    if exist('max_t','var') && (t1(i)>max_t || t1(i)<min_t)
        continue;
    end
    
    
    itmp=find(abs(t1(i)-t2)<=dtmax);
    
    itmp2=setdiff(itmp,i2);
    if ~isempty(itmp2)
        
        [val,itmp3]=min(abs(t1(i)-t2(itmp2)));
        i2(end+1)=itmp2(itmp3);
        i1(end+1)=i;
    end
 end

else
    
for i=1:length(t2)
    
    if exist('max_t','var') && (t2(i)>max_t || t2(i)<min_t)
        continue;
    end
    
    
    itmp=find(abs(t2(i)-t1)<=dtmax);
    
    itmp2=setdiff(itmp,i1);
    if ~isempty(itmp2) 
        
        [val,itmp3]=min(abs(t2(i)-t1(itmp2)));
        i1(end+1)=itmp2(itmp3);
        i2(end+1)=i;
    end
 end
end    
    
    
    

