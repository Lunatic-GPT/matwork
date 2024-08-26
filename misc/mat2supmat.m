function sF=mat2supmat(F,l,lr)
%sF=mat2supmat(F,l,lr)
%when lr = 'r',obtain the superoperator matrix for F as in v=F*u where u has a dimension of
%size(F,2)*l; Fsuper has dimension of size(F,1)*l by size(F,2)*l
% when lr = 'l', obtain the superoperator matrix for F as in v=u*F where u
% has a dimension of l*size(F,1); Fsuper has a dimension of size(F,2)*l by size(F,1)*l
% u and v are represented as column vectors u(:) and v(:) in the new space.

sz=size(F);

if strcmp(lr,'r')
sF=zeros(sz(1)*l,sz(2)*l);

for i=1:sz(1)
    for j=1:l
        row=i+(j-1)*sz(1);
        for alpha=1:sz(2)
         column=alpha+(j-1)*sz(2);
         sF(row,column)=F(i,alpha);
         
        end
    end
end

else
    
sF=zeros(sz(2)*l,sz(1)*l);

for i=1:sz(2)
    for j=1:l
        row=j+(i-1)*l;
        for alpha=1:sz(1)
         column=j+(alpha-1)*l;
         sF(row,column)=F(alpha,i);
         
        end
    end
end
end




