function  res = tTransform(m)

%res = tTransform(m)
% l2 norm of each row will be normalized to one.
% Implements an arbitrary transformation along time by m;
% m can be rank-deficient.
%
res.adjoint = 0;
m2=conj(m).*m;
m2=sqrt(sum(m2,2));
m2=repmat(m2,[1,size(m,2)]);

res.xfm = m./m2;

res = class(res,'tTransform');

