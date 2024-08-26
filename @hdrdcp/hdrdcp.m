function res = hdrdcp(ref,b)
% res = Wavelet(ref)
%
% implements a wavelet operator
%
% (c) Michael Lustig 2007

res.adjoint = 0;
ref=ref(:)-mean(ref);

n=length(ref);
xfm=[ones(n,1),ref];
m=eye(n);


for i=1:n
    beta=xfm\m(:,i);
    m(:,i)=m(:,i)-xfm*beta;
end

xfm2=orth(m);

b=diag(1./b(:));

xfm3=cat(2,xfm,xfm2)*b;

res.ixfm=inv(xfm3);

res.xfm=xfm3;
res = class(res,'hdrdcp');
