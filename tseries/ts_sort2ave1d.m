function [ts_mn,ts_std]=ts_sort2ave1d(ts,nTR_trial,ind_mean,ind_bs,exclude)
% [ts_mn,ts_std]=ts_sort2ave1d(ts,nTR_trial[,ind_mean,ind_bs,exclude])
%ind_mean: 0 return the averaged time course. > 0, then calculate the
%          average over the ind_mean time points. 1 based.

ts_tmp = reshape(ts,nTR_trial,length(ts)/nTR_trial);

if exist('exclude','var') && ~isempty(exclude)
ts_tmp(:,exclude) = [];
end

if exist('ind_bs','var') && ~isempty(ind_bs)
 ts_tmp = ts_tmp - repmat(mean(ts_tmp(ind_bs,:),1),[nTR_trial,1]);
end

if ~exist('ind_mean','var') || ind_mean(1) == 0 
 ts_mn=mean(ts_tmp,2); 
 ts_std = std(ts_tmp,0,2)/sqrt(size(ts_tmp,2));
else
 ts_tmp = mean(ts_tmp(ind_mean,:),1);
 ts_mn = mean(ts_tmp,2);
 ts_std = std(ts_tmp,0,2)/sqrt(size(ts_tmp,2));
end
    
