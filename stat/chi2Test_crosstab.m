function p = chi2Test_crosstab(y)
%p = chi2Test_crosstab(y)
%y is a crosstab

N=sum(y(:));

x1=zeros(N,1);
x2=zeros(N,1);

n=1;
for i=1:size(y,1)
    x1(n:n+sum(y(i,:))-1)=i;
    m=1;
    for j=1:size(y,2)
      x2(n-1+m:n-1+m+y(i,j)-1)=j;
      m=m+y(i,j);
    end
    n=n+sum(y(i,:));
end


[tab,chi2,p]=crosstab(x1,x2);





  