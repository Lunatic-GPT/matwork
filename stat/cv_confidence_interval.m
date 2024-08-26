function [res,res_relative]=cv_confidence_interval(cv_true,n,ci)
% ci: confidence interval; 
% sd: measured sd
%McKay, A.T., Distribution of the Coefficient of Variation and the Extended "t" Distribution. 
%Journal of the Royal Statistical Society, 1932. 95(4): p. 695-698.
if ~exist('ci','var')
    ci=0.95;
end
%cdf_cv(v, cv_true,n)
%{
func=@(x) cdf_cv(x,cv_true,n)-(1-ci)/2;
u=fzero(func,cv_true/2);
func=@(x) cdf_cv(x,cv_true,n)-(1+ci)/2;
l=fzero(func,sd*2);
res=[l,u];
%}
func=@(x) cdf_cv(cv_true+x,cv_true,n)-cdf_cv(cv_true-x,cv_true,n)-ci;

x=fzero(func,cv_true/2);

res=[cv_true-x,cv_true+x];

res_relative=res/cv_true;






