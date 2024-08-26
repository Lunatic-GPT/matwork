function res=scale2n(d,n,lmt)
% scale data to 1:n;

d=double(d);
if ~exist('lmt','var')
    lmt=[min(d(:)),max(d(:))];
end

r=lmt(2)-lmt(1);
step=r/n;

res=ceil((d-lmt(1))/step);

res(res<1)=1;
res(res>n)=n;

%res=uint16(res);