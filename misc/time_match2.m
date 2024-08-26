function [d1_new,d2_new]=time_match2(t1,t2,d1,d2,dtmax,max_t)
% find the times in t1 and t2 that are within dtmax of each other and
% return their indices
% if a time matches with more than one time in the other array, then return
% the closest time in the other array.
% assume t1 and t2 are both sorted in the same order.
d1_new=[];
d2_new=[];


for i=1:length(t1)
    
    if exist('max_t','var') && t1(i)>max_t
        continue;
    end
    
    ind = find(t1(i)<t2(2:end) & t1(i)>t2(1:end-1));
    
    if ~isempty(ind)
        tmp=d2(ind,:)*(t1(i)-t2(ind))+d2(ind+1,:)*(t2(ind+1)-t1(i));
        tmp=tmp/(t2(ind+1)-t2(ind));
        d2_new(end+1,:)=tmp;
        d1_new(end+1,:)=d1(i,:);
    elseif t1(i)>=t2(end) && t1(i)-t2(end)<=dtmax
            d2_new(end+1,:)=d2(end,:);
            d1_new(end+1,:)=d1(i,:);
    elseif t1(i)<=t2(1) && -t1(i)+t2(1)<=dtmax
            d2_new(end+1,:)=d2(1,:);
            d1_new(end+1,:)=d1(i,:);
    else
        
    end
    
end

        
    
    
    

