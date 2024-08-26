function res=roi_boundary(roi_all,nlayer)

unq=unique(roi_all(:));
unq(unq==0)=[];
res=0*roi_all;

for i=1:length(unq)
   
    roi=roi_all==unq(i);
    roi2=roi*0;
    
    for j=1:nlayer
        roi2=roi2|bwmorph(roi>0,'remove');
        
        roi=roi&roi2==0;
    end
    
    res(roi2>0)=unq(i);
    
end