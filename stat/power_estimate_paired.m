function pwr = power_estimate_paired(pval,diff,n1)
%  pwr = power_estimate_paired(pval,diff,n1)
% n1: number of pairs
% diff is the effect size, i.e. true difference normalized by SD of the difference (paired).
% same as pwr.t.test(type='one.sample')

    sd=1;
    try
        func=@(x) tcdf(x,n1-1)-pval/2;
        thr= abs(fzero(func,-1.96));  % normalized threshold
    catch
        disp([sd,n1,pval]);
    end
    
    diff=abs(diff); % sign does not matter; assume position change
  
    pwr=1-nctcdf(thr,n1-1,diff*sqrt(n1));
    