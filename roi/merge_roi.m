function roi=merge_roi(fname)
roi=load(fname);
roi=sum(roi.roi,4)>0;

save(fname,'roi');