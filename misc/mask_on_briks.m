function d=mask_on_briks(d,m,briks,val)

% d=mask_on_briks(d,m,briks,val)
% set briks of d with m>0 to value val;
% ndims(d)=ndims(m)+1
% briks selects the components along the last dimension.

sz=size(d);

d=reshape(d,[prod(sz(1:end-1)),sz(end)]);


tmp=d(:,briks);

tmp(repmat(m(:)>0,[1,length(briks)]))=val;

d(:,briks)=tmp;

d=reshape(d,sz);

