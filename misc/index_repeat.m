function ind=index_repeat(a)

n=0*a;
for i=1:length(a)
  n(i)=sum(a(i)==a);
end
ind=find(n>1);

    