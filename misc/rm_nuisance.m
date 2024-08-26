function ts = rm_nuisance(ts,nu,tau)
% ts = rm_nuisance(ts,nu)
%nu = [nu,ones(size(nu,1),1)];
b = zeros(1,2*tau+1);
for i=1:2*tau+1
    i1 = i-1-tau;
    
    if i1<1
        nu2 = nu(-i1+1:end,:);  
        ts2 = ts(1:end+i1);
        
    else 
         nu2 = nu(1:end-i1,:);
        ts2 = ts(i1+1:end);
    end
    
     b(i) = inv(nu2'*nu2)*nu2'*ts2;
end

[bmax,d] = max(b);

i1 = d-1-tau;
    
    if i1<1
        nu2 = nu(-i1+1:end,:);  
        ts(1:end+i1) = ts(1:end+i1) - nu2*bmax;
        
    else 
         nu2 = nu(1:end-i1,:);
         ts(i1+1:end) = ts(i1+1:end) - nu2*bmax;
    end
