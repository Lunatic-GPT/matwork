function res=getv_ind(d,ind)
% ind is the index for d(:,:,...,1): the first dimensioins should match.

sz=size(d);

d=reshape(d,[prod(sz(1:end-1)),sz(end)]);

res=d(ind,:);


