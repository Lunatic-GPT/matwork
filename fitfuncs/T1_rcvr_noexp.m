function y=T1_rcvr_noexp(b,x)
% y=b1*(1-b2*exp(-(x/b3)^b4));
y=b(1)*(1-b(2)*exp(-(x/b(3))^b(4)));

