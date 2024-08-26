function draw_color_bar(cm)
figure;
nc=size(cm,1);
z=repmat((1:nc)',[1,round(0.1*nc)]);


imshow(z,cm);

