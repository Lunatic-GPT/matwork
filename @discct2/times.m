function res = times(a,b)

if isa(a,'discct2') == 0
    error('In  A.*B only A can be discct2 operator');
end

    
    n=a.n;
    m=a.m;
   hx=zeros(size(b));
   for i=1:size(b,1)/m
     for j=1:size(b,2)/n
       if a.adjoint  
        hx((i-1)*m+1:(i-1)*m+m,(j-1)*n+1:(j-1)*n+n)=idct2(b((i-1)*m+1:(i-1)*m+m,(j-1)*n+1:(j-1)*n+n));
       else
           hx((i-1)*m+1:(i-1)*m+m,(j-1)*n+1:(j-1)*n+n)=dct2(b((i-1)*m+1:(i-1)*m+m,(j-1)*n+1:(j-1)*n+n));
       end
     end
   end

    res = hx;


