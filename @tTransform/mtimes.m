function res = mtimes(a,b)

sz=size(b);
b2=reshape(b,[prod(sz(1:end-1)),sz(end)]);
b2=permute(b2,[2,1]);

if a.adjoint
    
    res = a.xfm'*b2;
  
else

    res = a.xfm*b2;
  
end

res=permute(res,[2,1]);

sz2=sz;
sz2(end)=size(res,2);
res=reshape(res,sz2);





    
