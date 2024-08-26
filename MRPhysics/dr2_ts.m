function dr2=dr2_ts(d,te,bl,order)
%dr2=dr2_ts(d,te,bl[,order=0])
% order: 1 for first order baseline correction.
if ~exist('order','var')
    order=0;
end

s0=mean(d(bl));

if order==1
  d=reshape(d,[1,1,1,length(d)]);
  
  d=ts_detrend(d,bl,1);
  d=squeeze(d)+s0;
  d=d';
end    

dr2=-log(d/s0)/te;