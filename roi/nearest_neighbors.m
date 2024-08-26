function m2=nearest_neighbors(m)
%m=neighbors_mask(pos,dist,sz)
%pos are one-based integers; n*3
% dist 

ind=find(m>0);

m2=zeros(size(m));
 pos=ind2subb(size(m),ind); 
 
for i=1:length(ind)
  
    cb=-1:1;
 m2(pos(i,1)+cb,pos(i,2)+cb,pos(i,3)+cb)=1;
    
end

m2=m2&~m;