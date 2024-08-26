function m=Dymat(dim)
%m=Dymat(dim)
% a matrix to return the difference between neighboring elements along the columns
% dim is the number of rows
% y=y*m

m=zeros(dim,dim);

for i=2:dim
    m(i,i)=1;
    m(i-1,i)=-1;
end
m(1,1)=1;
m(end,1)=-1;



