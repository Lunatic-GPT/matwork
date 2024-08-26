function psf_cs(xfm,mask,ind)
%psf_cs(xfm,mask,ind)

sz=size(mask);
img=zeros(size(mask));
FT = p2DFT_rect(mask, sz, 1, 2);

if length(ind)==2
    ind(3)=1;
end

img(ind(1),ind(2),ind(3))=1;

img2=xfm'*img;
img3=FT*img2;
img4=FT'*img3;
img5=xfm*img4;

img_tmp=img5;

img_tmp(ind(1),ind(2),ind(3))=0;
mx=max(abs(img_tmp(:)));
disp(mx/abs(img5(ind(1),ind(2),ind(3))));

img6=abs(img5)/mx*100;
figure;mesh(img6,ones(sz));colormap(gray(100));
xlim([0,sz(1)]);
ylim([0,sz(2)]);

