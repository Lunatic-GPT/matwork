function [id,id_small]=index_within_limit(id,sz)
% id: n*2;
% sz: n*1 or 1*n
% the returned indices in id(i) will be between [1,sz(i)]
% id_small: n*2; the corresponding indices between 1 and
% (id(:,2)-id(:,1)+1)

n=size(id,1);
id_small=[ones(n,1),(id(:,2)-id(:,1)+1)];

for i=1:length(sz)
   if id(i,1)<1
       
       id_small(i,1)=-id(i,1)+2;
       id(i,1)=1;
       
   end
   
   if id(i,2)>sz(i)
       id_small(i,2)=id_small(i,2)-id(i,2)+sz(i);
        
       id(i,2)=sz(i);
       
   end
    
end