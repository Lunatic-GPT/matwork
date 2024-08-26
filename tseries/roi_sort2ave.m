function [mn,s] = roi_sort2ave(data,nTR_trial,mask,ts_ind,bs_ind)
% [mn,s] = roi_sort2ave(data,nTR_trial,mask,ts_ind[,bs_ind])
% mn and s are vectors containing mean and standard deviations of the
% voxels in mask>0
% ts_ind: 1 based. time range for average.
% bs_ind: 1 base. baseline
if isa(mask,'char')
    mask = BrikLoad(mask);
end

if isa(data,'char')
   ts = BrikLoad(data);
else 
    ts = data;
end

sz = size(ts);
ts = reshape(ts,[sz(1:3),nTR_trial,sz(4)/nTR_trial]);

if exist('bs_ind','var') && ~isempty(bs_ind)
    a = mean(ts(:,:,:,ts_ind,:),4)-mean(ts(:,:,:,bs_ind,:),4);
else 
    a = mean(ts(:,:,:,ts_ind,:),4);
end

a = squeeze(a);

mn = mean(a,4);
s = std(a,0,4)/sqrt(size(a,4));

mn = mn(mask>0);
s = s(mask>0);

