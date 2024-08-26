function [W,p] = kkc(x,m,noplot)
% W= kkc(x,m)
% x is a four dimensional data. m is a three dimensional mask
% calculate the Kendall's coefficient of concordance


sz = size(x);

nv = length(find(m>0));

m = repmat(m>0,[1,1,1,sz(4)]);

x = x(m);

ts = reshape(x,[nv,sz(4)]);
ts_rank = zeros(nv,sz(4));

for i=1:nv
    [ts_sort,ind]=sort(ts(i,:));
    ts_rank(i,ind) = 1:sz(4);
end

a = ts_rank';

%a_mn = mean(a,2);
%y = (a - repmat(a_mn,[1,nv])).^2;
%K = nv^2*sz(4)*(sz(4)^2-1)/12/(nv-1);
%W = 1 - 1/K*sum(y(:));

a_sum = sum(a,2);
S = sum((a_sum - mean(a_sum)).^2);
W = 12*S/nv^2/(sz(4)^3-sz(4));

p = chi2Test(W*(sz(4)-1)*nv,sz(4)-1);

if exist('noplot','var') && ~noplot
 figure;
 
   imagesc(ts_rank);
   colormap(gray);
  
end















