function ts_s=Tsmooth(ts,method)
% ts_s=Tsmooth(ts,method)
% if both dimensions are non-singlton, then the first dimension should be
% time.
if ndims(ts)>2
    error('Number of dimensions should be equal to 2');
end

ts2=ts;
if size(ts,1)==1
    ts2=ts';    
end
    ts2= cat(1,ts2(1,:),ts2,ts2(end,:));

    
if method ==1
   ts_s = 0.15*ts2(1:end-2,:)+0.7*ts2(2:end-1,:)+0.15*ts2(3:end,:);
elseif method ==2
   
    ts3= cat(3,ts2(1:end-2,:),ts2(2:end-1,:),ts2(3:end,:));
    a=min(ts3,[],3);
    b=max(ts3,[],3);
    c=median(ts3,3);
    ts_s = 0.15*a+0.7*c+0.15*b;  
else
    error('Unknown method');
end

if size(ts,1)==1
    ts_s=ts_s';    
end



