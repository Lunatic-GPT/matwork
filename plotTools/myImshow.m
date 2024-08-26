function myImshow(d,range,cm)


d=scale2n(d,size(cm,1),range);

imshow(d,cm);