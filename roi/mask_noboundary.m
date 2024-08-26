function m=mask_noboundary(m)

m1e=m([1,end],:);
m1e(m1e>0)=0;
m([1,end],:)=m1e;

m1e=m(:,[1,end]);
m1e(m1e>0)=0;
m(:,[1,end])=m1e;

