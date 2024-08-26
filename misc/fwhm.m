function res = fwhm(sp)

[pk,ind]=max(sp);


f1=@(x) interp1(1:ind,sp(1:ind),x,'linear',pk)-pk/2;
x1=fzero(f1,ind-5);

f2=@(x) interp1(ind:length(sp),sp(ind:end),x,'linear',pk)-pk/2;
x2=fzero(f2,ind+10);

res=x2-x1;
