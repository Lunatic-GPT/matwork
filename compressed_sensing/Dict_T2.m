function D=Dict_T2(K,n,te,R2)

% K: sparsity
% n: number of atoms
% te: te values
% T2 values in the dictionary

nte=length(te);
nR2=length(R2);


x=zeros(nte,nR2);

for i=1:nte
    for j=1:nR2
           x(i,j)=exp(-te(i)*R2(j));  
    end
end


y = x*diag(1./sqrt(sum(x.*x)));


params.data = y;
params.Tdata = K;
params.dictsize = n;
params.iternum = 30;
params.memusage = 'high';

D = ksvd(params,'');





