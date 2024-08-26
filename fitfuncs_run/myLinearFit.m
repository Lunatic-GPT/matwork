function [k,b]=myLinearFit(d,x)
% size(d,ndim)=length(x)

sz=size(d);
d=reshape(d,[prod(sz(1:end-1)),sz(end)]);

xx=[x(:),ones(length(x),1)];

res=xx\permute(d,[2,1]);

res=permute(res,[2,1]);

k=reshape(res(:,1),sz(1:end-1));
b=reshape(res(:,2),sz(1:end-1));



