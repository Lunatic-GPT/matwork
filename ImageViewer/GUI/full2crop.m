function img=full2crop(f,rect_pos)

img=f(rect_pos(2):(rect_pos(2)+rect_pos(4)),rect_pos(1):(rect_pos(1)+rect_pos(3)));

