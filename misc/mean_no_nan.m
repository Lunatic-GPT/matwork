function res=mean_no_nan(d,dim)

sel=isnan(d);

d(isnan(d))=0;

res=mean(d,dim)*size(d,dim)/(size(d,dim)-sum(sel,dim));



