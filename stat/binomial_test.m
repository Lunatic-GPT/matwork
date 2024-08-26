function p=binomial_test(n,N, prob)
% prob is a value between 0 and 1; the prob of observing case A
% returns the p value of observing <=n case A out of N observations

p=0;
for i=0:n
   p=p+CnN(i,N)*prob^n*(1-prob)^(N-n); 
end



function  res=CnN(n,N)
res=factorial(N)/factorial(n)/factorial(N-n);

