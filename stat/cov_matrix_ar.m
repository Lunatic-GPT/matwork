function res=cov_matrix_ar(n,a)
%res=cov_matrix_ar(n,a)
res=zeros(n,n);

for i=1:n
    for j=1:n
        res(i,j)=a^abs((i-j));
    end
end