function res=sampleSize_SD_estimate(ci,pe,n0)
% ci: confidence interval; 
% pe: percent error
% reference: Greenwood and Sandomire, Journal of American Statistical
% Association, 45: 257-260.

func=@(x) ci_func(pe,x)-ci;
res=fzero(func,n0);

function res=ci_func(pe,n)
p1=chi2cdf(n*(1+pe)^2,n);
p2=chi2cdf(n*(1-pe)^2,n);

res=p1-p2;


