function m=Dxmat(dim)
%m=Dxmat(dim)
% a matrix to return the difference between neighboring elements along the rows
% dim is the number of rows
% dx=m*x
m=zeros(dim,dim);

for i=2:dim
    m(i,i)=1;
    m(i,i-1)=-1;
end
m(1,1)=1;
m(1,end)=-1;


