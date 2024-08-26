function PVS_orientation_colormap(d,nrow,ncol,roi)
% 


d=ri(d);
d=squeeze(abs(d));
 
%clr_roi=lines(4);

clr_roi=ones(4,3)*255;
 
d=flip(permute(d,[2,1,3,4]),1);

roi=flip(permute(roi,[2,1,3,4]),1);

d=d/max(d(:));
d=uint8(d*2024);
d=add_roi2imageRGB(d,roi,clr_roi);

sz=size(d);

im=reshape(d,[sz(1),sz(2),ncol,nrow,sz(4)]);
im=permute(im,[1,4,2,3,5]);
im=reshape(im,[sz(1)*nrow,sz(2)*ncol,sz(4)]);



figure;imshow(im);



