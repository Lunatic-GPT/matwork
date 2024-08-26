function roi2=roi_subset_exc(roi,iexc)
% exclude pixels with values in ikeep
roi2=roi;
for i=1:length(iexc)
   roi2(roi==iexc(i))=0;
   
end