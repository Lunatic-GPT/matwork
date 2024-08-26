function res=sampleSize_estimate(pwr,pval,diff,ratio)

% ratio: the number of subjects between group 1 and group 2: default 1;
% 

if ~exist('ratio','var')
    ratio=1;
end

func=@(n) (power_estimate_paired(pval,diff,ratio)-pwr).^2;

res=fminsearch(func,20);

