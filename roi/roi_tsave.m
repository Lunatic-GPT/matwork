function o = roi_tsave(data,mask,ts_ind,bs_ind)
%o = roi_tsave(data,mask,ts_ind[,bs_ind])
% ts_ind,bs_ind: both 1 based.
% if ts_ind = 0, then return the whole time course in the second dimension;
% first dimension is voxel.

if isa(mask,'char')
    mask = BrikLoad(mask);
end

if isa(data,'char')
   ts = BrikLoad(data);
else 
    ts = data;
end

if exist('bs_ind','var') && ~isempty(bs_ind)
  ts = ts - repmat(mean(ts(:,:,:,bs_ind),4),[1,1,1,size(ts,4)]);
end

if ts_ind > 0
  dp = mean(ts(:,:,:,ts_ind),4);
  o = dp(mask>0);
else
    
  nv = length(find(mask>0));  
  nm = repmat(mask,[1,1,1,size(ts,4)]);
  ts2d = ts(nm>0);
  o = reshape(ts2d,[nv,size(ts,4)]);
  
end 
