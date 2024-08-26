function  res = KLT(data)
% KLT(data)

res.adjoint = 0;

v=KLT_matrix(data);
%[v,d]=eig(cov);

res.m=v';
res = class(res,'KLT');



