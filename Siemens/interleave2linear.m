function d2=interleave2linear(d,dim)
% for odd number of slices, the output order is 1,(n+1)/2+1,2,(n+1)/2+2,...(n-1)/2,(n+1)/2
% for even number of slices, the output order is n/2+1,1,n/2+2,2,...n,n/2
% for axial slices, the first slice is on the feet side

ns=size(d,dim);

ndim=ndims(d);

arr=1:ndim;

arr([3,dim])=arr([dim,3]);

d=permute(d,arr);

d2=d;

if mod(ns,2)==0
    d2(:,:,2:2:end,:,:)=d(:,:,1:ns/2,:,:);
    d2(:,:,1:2:end-1,:,:)=d(:,:,ns/2+1:end,:,:);

else
    
    d2(:,:,1:2:end,:,:)=d(:,:,1:(ns+1)/2,:,:);
    d2(:,:,2:2:end-1,:,:)=d(:,:,(ns+1)/2+1:end,:,:);

end

    d2=permute(d2,arr);

    
    