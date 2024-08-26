function ind=sub2indb(mtrx,sub)

ind=zeros(1,size(sub,1));
for j=1:size(sub,1)
ind(j)=sub(j,1);
for i=2:size(sub,2)
    ind(j)=ind(j)+(sub(j,i)-1)*prod(mtrx(1:i-1));
end

end