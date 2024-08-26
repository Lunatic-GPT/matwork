function p = normalityTest(d)

% fdist_nc(x,n,m,ld)

mn=mean(d);

sd=std(d);
x=(d-mn)/sd;
[h,p]=kstest(x,'Alpha',0.01);


figure;
cdfplot(x);
hold on
x_values = linspace(min(x),max(x));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best');


