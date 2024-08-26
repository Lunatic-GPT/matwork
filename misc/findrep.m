function res=findrep(a)

[v,ind]=unique(a);

ind2=setdiff(1:length(a),ind);
res=a(ind2);