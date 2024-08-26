function d = fdist_nc(x,v1,v2,ld)
% fdist_nc(x,n,m,ld)

d = zeros(size(x));


for i=1:length(d)
    k = 0;
while 1
    if x(i) <=0
        break;
    end
      
    tmp = exp(-ld/2)*(ld/2)^k/beta(v2/2,v1/2+k)/gamma(k+1);
    tmp2 = (v1/v2)^(v1/2+k);
    tmp3 = (v2/(v2+v1*x(i)))^((v1+v2)/2+k)*x(i)^(v1/2-1+k);
    
    inc = tmp*tmp2*tmp3;
    if inc < d(i)/1000
        break;
    end
      d(i) = d(i) + inc;
      k=k+1;
end
end