function d2=interleave2linear_times(ns,mb_factor,tr)
% for odd number of slices, the output order is 1,(n+1)/2+1,2,(n+1)/2+2,...(n-1)/2,(n+1)/2
% for even number of slices, the output order is n/2+1,1,n/2+2,2,...n,n/2
% for axial slices, the first slice is on the feet side


ns2=ns/mb_factor;

torder=zeros(ns2,1);

if mod(ns2,2)==0
    torder(2:2:end)=1:ns2/2;
    torder(1:2:end-1)=ns2/2+1:ns2;
else
    torder(1:2:end)=1:(ns2+1)/2;
    torder(2:2:end-1)=(ns2+1)/2+1:ns2;
end

torder=repmat(torder,[1,mb_factor])-1;

disp([(1:ns)',round(torder(:))+1]);

%disp([(1:ns)',round(torder(:)*tr/ns2)]);





    
    