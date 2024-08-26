function res=std_no_nan(d,flag,dim)

sel=isnan(d);

d(isnan(d))=0;

res=std(d,flag,dim)*size(d,dim)/(size(d,dim)-sum(sel,dim));



