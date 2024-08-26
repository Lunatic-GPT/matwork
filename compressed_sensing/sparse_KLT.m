function sparse_KLT(d,t)


N=length(t);
sigma=zeros(N,N);
for i=1:N
 
    
 t=t(:);

[V,D]=eig(cov(t));
KLT=V'*t';