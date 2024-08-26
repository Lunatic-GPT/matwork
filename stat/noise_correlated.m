function x=noise_correlated(cov,N)

[v,d]=eig(cov);

nr=randn_white([size(v,1),N]);

ni=randn_white([size(v,1),N]);

%x=real(v)*sqrt(d)*nr+imag(v)*sqrt(d)*ni*1i;

%x=v*sqrt(d)*nr;
x=v*sqrt(d)*(nr+1i*ni);
x=conj(x);