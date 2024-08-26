function PValue = FisherExactTest_matrix(y)
% p = chi2Test(chi2,dof)
% the value from this function seems wrong. Diff from fisher.test in R and
% fishertest in Matlab
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


[ Sig,PValue,ContigenMatrix] =FisherExactTest(x1,x2);





  