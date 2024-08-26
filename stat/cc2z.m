function z=cc2z(cc,n)
% [p,t]=cc2tp(cc,n)
%cc, cc value. n. the number of time points in the time course to calculate
%the time course

t = sqrt(n-2)*cc/sqrt(1-cc^2);

z=t2z(t,n-2);

disp(0.5*(log(1+cc)-log(1-cc)));
p = tTest(n-2,t);
%p = (1-cdf('t',t,n-2))*2;

