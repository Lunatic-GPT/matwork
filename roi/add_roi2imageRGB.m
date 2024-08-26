function img=add_roi2imageRGB(img,rois,clr_roi)
%img is 4D with size(img,4)=3
%rois is 3D with different values for different ROI;
%clr_roi: n*3;

val=unique(rois(:));
val(val==0)=[];

for i=1:size(img,3)

    for j=1:length(val)
       
       roic=bwmorph(rois(:,:,i)==val(j),'remove');

       img(:,:,i,1)=setv_roi(img(:,:,i,1),roic,clr_roi(val(j),1));
       img(:,:,i,2)=setv_roi(img(:,:,i,2),roic,clr_roi(val(j),2));
       img(:,:,i,3)=setv_roi(img(:,:,i,3),roic,clr_roi(val(j),3));
        
    end


end



