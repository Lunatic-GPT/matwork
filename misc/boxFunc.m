function res=boxFunc(l,u,x)

res=ones(size(x));
res(x<l | x>u)=0;
