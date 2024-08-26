function [img,cm]=draw_image_roi(img,scale,rois)

% [img,cm]=draw_image_roi(img,scale,varargin)
% rois give the different rois that will be ploted in different colors

img=double(img);
if isempty(scale)
scale=min_max(img(:));
end


img=(img-scale(1))/diff(scale)*100;

img(img>100)=100;
cm=gray(100);
clr=[1,0,0;0,1,0;0,0,1;1,1,0;1,0,1];

for i=1:length(rois)
    
roic=bwmorph(rois{i}>0,'remove');
[img,cm]=combine_over_under(img,roic,cm,clr(i,:),roic);

end

img=scaled2rgb(img,cm);

if nargout==0

imshow(img);

end



