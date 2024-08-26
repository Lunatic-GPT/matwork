function [mn,sem] = groupAverage_rdk(data,coh4)
%[mn,sem] = groupAverage_rdk(data,coh4)
% data is assumed to have dimensions of [npts,4,narea,nsubjects];
% or [4,narea,nsubjects];
sz = size(data);

mn = zeros(sz(1:end-1));
sem = zeros(sz(1:end-1));
if length(sz) ==4 && sz(2) == 4
    
    mn(:,1,:) = mean(data(:,1,:,coh4),4);
    sem(:,1,:) = std(data(:,1,:,coh4),0,4)/sqrt(sz(end));
    mn(:,2:4,:) = mean(data(:,2:4,:,:),4);
    sem(:,2:4,:) = std(data(:,2:4,:,:),0,4)/sqrt(sz(end));
    
elseif length(sz) == 3 && sz(1) == 4
    
    mn(1,:) = mean(data(1,:,coh4),3);
    sem(1,:) = std(data(1,:,coh4),0,3)/sqrt(sz(end));
    mn(2:4,:) = mean(data(2:4,:,:),3);
    sem(2:4,:) = std(data(2:4,:,:),0,3)/sqrt(sz(end));
   
else
    error('dimension error');
end

   
    




 