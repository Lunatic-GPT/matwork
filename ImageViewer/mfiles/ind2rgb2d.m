function rgb=ind2rgb2d(img,cm)
rgb = zeros([size(img),3]);
for i=1:size(img,1)
    for j=1:size(img,2)
        rgb(i,j,:)=ind2rgb(round(img(i,j)),cm);
    end
end
        