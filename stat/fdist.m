function d=fdist(x,n,m)
% d=fdist(x,n,m)
% calculate the probability density of a f distribution

d = zeros(size(x));

for i=1:length(x)
%  d(i) = gamma((n+m)/2)*n^(n/2)*m^(m/2)/gamma(n/2)/gamma(m/2);

%  d(i) = d(i)*x(i)^(n/2-1)/(m+n*x(i))^((n+m)/2);

tmp = (n*x(i))^(n-2)*m^m/(n*x(i)+m)^(m+n);
tmp = sqrt(tmp);
d(i) = tmp/beta(n/2,m/2);
end

