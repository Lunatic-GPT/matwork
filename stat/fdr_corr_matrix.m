function pmask=fdr_corr_matrix(p,q)

pmask=0*p;

m=p*0;
for i=1:size(p,1)
    m(i,1:i-1)=1;
end

ind=find(m(:));
p2=p(ind);
[tmp,i_ind]=fdr(p2,q);

pmask=p<=tmp;



