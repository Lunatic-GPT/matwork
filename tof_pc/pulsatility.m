function [res,eres]=pulsatility(ts,ind_min,ind_max,method)
% ts is npt*nvessels
% ind_min and ind_max are time point indices for calculating pulsatility
% method 1: average before calculate pulsatility (default)
% 2: calculat pulsatility before average

if ~exist('method','var')
    method=1;
end

if method==1
   
 res= mean(mean(ts(ind_max,:),1)-mean(ts(ind_min,:),1))./mean(mean(ts(ind_max,:),1)+mean(ts(ind_min,:),1))*2;

eres=NaN;
    
else
p= (mean(ts(ind_max,:),1)-mean(ts(ind_min,:),1))./(mean(ts(ind_max,:),1)+mean(ts(ind_min,:),1))*2;

res=mean(p);
eres=std(p)/sqrt(size(p,2));


end

