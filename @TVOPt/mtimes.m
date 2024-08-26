function res = mtimes(a,b)

if size(b,4)==1
        res=0;
        return;
end
    
if a.adjoint
      
res = b(:,:,:,[1,1:end-1]) - b;
res(:,:,:,1) = -b(:,:,:,1);
res(:,:,:,end) = b(:,:,:,end-1);

else
 
res = b(:,:,:,[2:end,end]) - b;
end




    
