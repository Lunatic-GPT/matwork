function [mn,sem] = cmean_std(data,dim,mode,flag)
% [mn,sem] = cmean(data,dim[,mode=0])
% average over data along dimension dim.
% excluding data points which are equal to -9999.
% mode 0: only data points equal to -9999 are censored
% mode 1: exclude all data points with the same index along dimension
% dim.
% return NaN is no data points left after censoring.
% flag: flag in std calculated. Same as in matlab std function.
if ~exist('mode','var')
    mode = 0;
end

if ~exist('flag','var')
    flag = 0;
end

sz = size(data);
pd = 1:ndims(data);
pd(dim) =[];
pd = [dim,pd];
d = permute(data,pd);
d2d = reshape2d(d);

   mn = zeros(size(d2d,2),1);
   sem = zeros(size(d2d,2),1);
if mode ==0
   
  for i=1:size(d2d,2)
      
  b = d2d(d2d(:,i)>-9998,i);
  mn(i)=mean(b);
  sem(i) = std(b,flag)/sqrt(length(b));
  end
else
  incl = ~any(d2d<-9998,2);
  mn = mean(d2d(incl,:),1);
  sem = std(d2d(incl,:),flag,1)/sqrt(length(find(incl)));
end

if length(pd)>2
    sz(dim) =[];
    mn = reshape(mn,sz);
    sem = reshape(sem,sz);
else
    sz(dim) = 1;
    mn = reshape(mn,sz);
    sem = reshape(sem,sz);
end


