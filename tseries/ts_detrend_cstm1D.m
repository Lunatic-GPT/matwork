function ts=ts_detrend_cstm1D(ts,bs,order)
% ts=ts_detrend1D(ts,bs,order)
% both ts and bs should be column vectors.
% Order: polynomial baseline in addition to bs.  Use -1 to use only bs as
% baseline.
if order>=0
    x = (1:length(ts))';
    for i=0:order
      bs = [bs,x.^i];
    end
end
b = inv(bs'*bs)*bs'*ts;
ts = ts - bs*b;