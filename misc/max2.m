function [res,ind]=max2(d)

[res,ind]=max(d(:));
ind=ind2subb(size(d),ind);
