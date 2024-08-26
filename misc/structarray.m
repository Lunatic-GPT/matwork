function res=structarray(a,member,catd)

if ~exist('catd','var')
    catd=1;
end

res=[];
for i=1:length(a(:))   
    tmp=getfield(a(i),member);
    res=cat(catd,res,tmp);
end

