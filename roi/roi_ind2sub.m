function ind=roi_ind2sub(roi)

ind=ind2subb(size(roi),find(roi(:)>0));