function [x,y]=isocontour2D(d,val)

d2=repmat(d,[1,1,10]);
v=isosurface(d2,val);

i1=find(v.vertices(:,3)==9);
x=v.vertices(i1,1);
y=v.vertices(i1,2);

[x_sort,ind]=sort(x);
y_sort=y(ind);


if ~any(y_sort(1:end-1)>y_sort(2:end))  ||  ~any(y_sort(1:end-1)<y_sort(2:end))
  x=x_sort;
  y=y_sort;
else
  phi=atan2(x-mean(x),y-mean(y));
  [tmp,ind]=sort(phi);
  x=x(ind);
  y=y(ind);
end


