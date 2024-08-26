function [m,cost1,cost2]=cs_recon2(y,us,phi,lambda)
% [m,cost1,cost2]=cs_recon2(y,ufftm,phi,lambda)
% 2-D compressed sensing.
% y: undersampled data
% us: undersampling matrix.  same size as y. if set to true, the
% corresponding value in y will be used.
% phi:   sparsifying transform
% lambda: tuning parameter

y(~us)=0;
TolGrad=1e-3;
MaxIter=5;
alpha=0.1;
beta=0.6;

k=0;
m=ifft2(y);
%m=zeros(size(y));
g=gradf(m,us,y,lambda,phi);
dm=-g;

while k<MaxIter && g(:)'*g(:)>TolGrad
   t=1;
   while 1
       c1=cs_cost(m+t*dm,us,y,phi,lambda);
       c2=cs_cost(m,us,y,phi,lambda)+alpha*t*real(g(:)'*dm(:));
       if c1<=c2
           break;
       end
       t=beta*t;    
   end
   
   m=m+t*dm;
   gnew=gradf(m,us,y,lambda,phi);
   gamma=(gnew(:)'*gnew(:))/(g(:)'*g(:));
   g=gnew;
   dm=-g+gamma*dm;
   k=k+1;
    disp([k,gamma,c1,c2]);
end

[c,cost1,cost2]=cs_cost(m,us,y,phi,lambda);

function [c,c1,c2]=cs_cost(m,us,y,phi,lambda)

tmp=(sparse_FT(m,us,false)-y);
c1=sum(tmp(:)'*tmp(:));
tmp=phi(m);
c2=sum(abs(tmp(:)));
c=c1+lambda*c2;



function out=gradf(m,us,y,lambda,phi)


pm=phi(m);
mu=1e-9;

w=sqrt(abs(pm).^2+mu);

pm=pm./w;
ft1=sparse_FT(m,us,false)-y;

out=2*sparse_FT(ft1,us,true)+lambda*phi(pm');

%out=phi(pm);

