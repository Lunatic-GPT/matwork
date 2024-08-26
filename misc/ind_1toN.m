function res=ind_1toN(ind)

ind2=unique(ind);

res=0*ind;
for i=1:length(ind2)
   
    res(ind==ind2(i))=i;
end