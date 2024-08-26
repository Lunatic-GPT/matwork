function [m,cost1,cost2]=cs_recon(y,ufftm,phi,lambda)
% [m,cost1,cost2]=cs_recon(y,ufftm,phi,lambda)
% y: undersampled data
% ufftm: undersampled fft matrix, to convert m to y
% phi:   sparsifying transform
% lambda: tuning parameter

TolGrad=1e-3;
MaxIter=40;
alpha=0.05;
beta=0.6;
n=size(ufftm,2);

k=0;
m=zeros(n,1);

g=gradf(ufftm,m,y,lambda,phi);
dm=-g;

while k<MaxIter && sum(g'*g)>TolGrad
   t=1;
   
   while 1
       c1=cs_cost(ufftm,m+t*dm,y,phi,lambda);
       c2=cs_cost(ufftm,m,y,phi,lambda)+alpha*t*real(g'*dm);
       if c1<c2
           break;
       end
       t=beta*t;    
   end
   
   m=m+t*dm;
   gnew=gradf(ufftm,m,y,lambda,phi);
   gamma=sum(gnew'*gnew)/sum(g'*g);
   g=gnew;
   dm=-g+gamma*dm;
   k=k+1;
    disp(k);
end

[c,cost1,cost2]=cs_cost(ufftm,m,y,phi,lambda);



function [c,c1,c2]=cs_cost(ufftm,m,y,phi,lambda)


tmp=(ufftm*m-y);
c1=sum(tmp'*tmp);
c2=sum(abs(phi*m));
c=c1+lambda*c2;



function out=gradf(ufftm,m,y,lambda,phi)


pm=phi*m;
mu=1e-9;

    w=diag(sqrt(abs(pm).^2+mu));
out=2*(ufftm)'*(ufftm*m-y)+lambda*phi'*inv(w)*phi*m;

