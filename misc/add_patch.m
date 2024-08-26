function add_patch(x,y)

xp = [x(1),x(1),x(2),x(2)];
yp = [y(1),y(2),y(2),y(1)];
h=patch(xp,yp,[0.5,0.5,0.5]);
alpha(h,0.5);