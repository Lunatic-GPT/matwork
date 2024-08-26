function res=mean_omitnan2(d,dim1,dim2)

% omit entries if they are NaN for any index along dim2


sel=any(isnan(d),dim2);

sz=ones(1,ndims(d));
sz(dim2)=size(d,dim2);

sel=repmat(sel,sz);

d(sel)=NaN;

res=mean(d,dim1,'omitnan');

