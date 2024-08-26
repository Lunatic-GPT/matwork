function y=T1_sat_rcvr(b,x)
% y=b1*(1-exp(-x/b2));
% %4.2f[1-exp(-x/%3.2f)]
y=b(1)*(1-exp(-x/b(2)));

