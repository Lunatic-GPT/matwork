function c=p2cc(p,n)
% Calcuate the cc value given the p value and the number of data points.
% c=p2cc(p,n)
t = icdf('t',1-p/2,n-2);
c = t/sqrt(n-2+t*t);