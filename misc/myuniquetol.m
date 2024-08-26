function res=myuniquetol(a,tol)

b=unique(a);
res=b(1);


for i=2:length(b)
    
    
    if any(abs(res-b(i))<b(i)*tol)
        continue;
    end
    
    res(end+1)=b(i);
end
