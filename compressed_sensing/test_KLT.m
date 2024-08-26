r=zeros(100,1);
r(1:10)=1;
r(21:30)=1;
r(41:50)=1;
r(61:70)=1;
r(81:90)=1;
a=r;
%a=r';
ma=mean(a);
N=100;
c=zeros(N,N);
for i=1:N
    for j=1:N
     c(i,j)=(a(i)-ma)*(a(j)-ma);
    end
end


res=KLT(c);
b=a-mean(a);
b=reshape(b,[1,1,1,length(b)]);
y=res*b;
%[v,d]=eig(c);

%y=v'*(a-mean(a));

figure;plot(squeeze(y));
