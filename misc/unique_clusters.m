function [res,roi_val]=unique_clusters(roi)
% roi values from 1 to n
% roi_val: corresponding values in the original roi

roi_val=unique(roi(:));

roi_val(roi_val==0)=[];

res=0*roi;
for i=1:length(roi_val)
    res(roi==roi_val(i))=i;
    
end