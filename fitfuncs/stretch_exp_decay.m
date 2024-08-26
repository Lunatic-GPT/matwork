function y=stretch_exp_decay(b,x)
% y=b1*exp(-(x/b2)^b3);
y=b(1)*exp(-(x/b(2)).^b(3));

