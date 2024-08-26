function res=cdf_cv(v, cv_true,n)
%McKay  https://www.jstor.org/stable/2342041
%n: the sample size;

a=sqrt(n)/cv_true;

fun=@(x) pdf_cv(n,x,a);
res=integral(fun,0,v);





function res=pdf_cv(n,v_vec,a)

res=0*v_vec;
for i=1:length(v_vec)
v=v_vec(i);
    opv2=(1+v.*v);
nom=v.^(n-2).*exp(-a^2*v.^2/2./opv2)*gamma(n);
denom=2^((n-3)/2)*gamma((n-1)/2)*opv2.^(n/2)*sqrt(2*pi);

int=@(x) x.^(n-1)/gamma(n).*(exp(-1/2*(x-a./sqrt(opv2)).^2)+exp(-1/2*(x+a./sqrt(opv2)).^2));

y=integral(int,0,inf);
res(i)=y*nom./denom;

end
