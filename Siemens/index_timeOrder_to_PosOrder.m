function d=index_timeOrder_to_PosOrder(i,ns)
% for odd number of slices, the output order is 1,(n+1)/2+1,2,(n+1)/2+2,...(n-1)/2,(n+1)/2
% for even number of slices, the output order is n/2+1,1,n/2+2,2,...n,n/2
% for axial slices, the first slice is on the feet side

ind=zeros(1,ns);
if mod(ns,2)==0
    
    ind(1:ns/2)=2:2:ns;
    ind(ns/2+1:end)=1:2:ns-1;
    
else
    
    ind(1:(ns+1)/2)=1:2:ns;
    ind((ns+1)/2+1:end)=2:2:ns-1;
    
end

d=ind(i);