function res = mtimes(a,b)

    sz=size(b);
    b=reshape(b,[prod(sz(1:end-1)),sz(end)]);

    res=zeros(size(b));
        
    for i=1:size(b,1)
                if a.adjoint
                   res(i,:) = a.m'*b(i,:)';
                else
                   res(i,:) = a.m*b(i,:)';
                end
    end
    
    res=reshape(res,sz);





    
