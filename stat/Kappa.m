function res=Kappa(a,b)
% a,b: categories - values in 1:n
if isa(a,'cell')
    ab=cell2int([a(:);b(:)]);
    ab=reshape(ab,[length(ab)/2,2]);
    a=ab(:,1);
    b=ab(:,2);
end

p0=sum(a==b)/length(a);

c=unique([a(:);b(:)]);

pe=0;
for i=1:length(c)
pe=pe+sum(a==c(i)) *sum(b==c(i))/length(a)^2;

end

res=(p0-pe)/(1-pe);

function res=cell2int(c)

c2=unique(c);

for i=1:length(c)
   res(i)= strmatch(c{i},c2,'exact');
end
    
