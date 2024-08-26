function y=T1_rcvr_abs(b,x)
% y=b1*abs(1-b2*exp(-x/b3));
% %4.2f[1-%3.2fexp(-x/%3.2f)]
y=b(1)*abs(1-b(2)*exp(-x/b(3)));

