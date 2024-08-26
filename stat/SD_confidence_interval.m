function res=SD_confidence_interval(sd,n,ci)
% ci: confidence interval; 
% sd: measured sd
% reference: Greenwood and Sandomire, Journal of American Statistical
% Association, 45: 257-260.
if ~exist('ci','var')
    ci=0.95;
end

func=@(x) p1_func(sd/x,n)-(1-ci)/2;
u=fzero(func,sd/2);

func=@(x) p1_func(sd/x,n)-(1+ci)/2;
l=fzero(func,sd*2);

res=[l,u];

%ratio
func=@(r) p1_func(1+r,n)-p1_func(1-r,n)-ci;
r=fzero(func,1/2);
res(2,:)=sd*[1-r,1+r];


function p1=p1_func(ratio,n)

p1=chi2cdf(n*ratio^2,n);




