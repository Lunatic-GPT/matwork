function [a,w,normr]=deconvolve(ts,irf,wlimit)
%[a,w,normr]=deconvolve(ts,irf,wlimit)
if size(ts,1) ==1
    ts = ts';
end

if size(irf,1) ==1
    irf = irf';
end

normr_vec = zeros(1,wlimit);
avec = zeros(1,wlimit);
nts = length(ts);
nirf = length(irf);

for i=1:wlimit
    ref = zeros(nts,1);
    ref(1:i+nirf-1)=conv(ones(i,1),irf);
    [p,s]=polyfit(ref(1:nts),ts,1); 
    avec(i) = p(1);
    normr_vec(i) = s.normr;
end

[normr,w]=min(normr_vec);
a = avec(w);
