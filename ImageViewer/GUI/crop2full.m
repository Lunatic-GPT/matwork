function img=crop2full(c,rect_pos,f)
img = zeros(size(f));

img(rect_pos(2):(rect_pos(2)+rect_pos(4)),rect_pos(1):(rect_pos(1)+rect_pos(3)))=c;




