function o = roi_tsave(data,mask,ts_ind,bs_ind)
%o = roi_tsave(data,mask,ts_ind[,bs_ind])

if isa(mask,'char')
    mask = BrikLoad(mask);
end

if isa(data,'char')
   ts = BrikLoad(data);
else 
    ts = data;
end

if exist('bs_ind','var') && ~isempty(bs_ind)
  dp = mean(ts(:,:,:,ts_ind),4) - mean(ts(:,:,:,bs_ind),4);
else 
  dp = mean(ts(:,:,:,ts_ind),4);
end
o = dp(mask>0);

