function res=kgrid2pos(sz)

n=length(sz);

sz(length(sz)+1:3)=1;

    
x=get_vec(sz(1));
y=get_vec(sz(2));
z=get_vec(sz(3));

[xx,yy,zz]=meshgrid2(x,y,z);

res=[xx(:),yy(:),zz(:)];

res=res(:,1:n);




function x=get_vec(sz)

if mod(sz,2)==0
    x=-sz(1)/2:sz(1)/2-1;
    x=x/sz(1);
else
    x=-(sz-1)/2:(sz-1)/2;
    x=x/sz(1);
end

x=x(:);



