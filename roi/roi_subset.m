function roi2=roi_subset(roi,ikeep)

roi2=0*roi;
for i=1:length(ikeep)
   roi2(roi==ikeep(i))=ikeep(i);
   
end