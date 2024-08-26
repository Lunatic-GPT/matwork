function [p,t]=cc2tp(cc,n)
% [p,t]=cc2tp(cc,n)
%cc, cc value. n. the number of time points in the time course to calculate
%the time course

t = sqrt(n-2)*cc/sqrt(1-cc^2);
p = tTest(n-2,t);
%p = (1-cdf('t',t,n-2))*2;

