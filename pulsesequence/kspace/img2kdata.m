function res = img2kdata(img,voxsize,k)


sz=size(img);

x=-sz(1)/2:sz(1)/2-1;
y=-sz(2)/2:sz(2)/2-1;
z=-sz(3)/2:sz(3)/2-1;

x=x*voxsize(1);
y=y*voxsize(2);
z=z*voxsize(3);


[x,y,z]=ndgrid(x,y,z);

res=img.*exp(-1i*(x*k(1)+y.*k(2)+z.*k(3)))*sinc(k(1)*voxsize(1)/pi)*sinc(k(2)*voxsize(2)/pi)*sinc(k(3)*voxsize(3)/pi)*prod(voxsize);


res=mean(res(:));