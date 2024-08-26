function [V,mat]=KLT_matrix(d)

nd=ndims(d);
sz=size(d);

d=reshape(d,[prod(sz(1:nd-1)),sz(end)]);

N=sz(end);
mat=zeros(N,N);

dm=mean(d,2);
%dm=zeros(size(d,1),1);
for i=1:N
    for j=i:N
     mat(i,j)=mean(conj(d(:,i)-dm).*(d(:,j)-dm));
     mat(j,i)=conj(mat(i,j));
    end
end

[V,D]=eig(mat);
V=fliplr(V);