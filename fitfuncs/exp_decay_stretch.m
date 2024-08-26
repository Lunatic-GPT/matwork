function y=exp_decay_stretch(b,x)
% y=b1*exp(-x/b2);
% %4.3f*exp(-x/%4.3f)
y=b(1)*exp(-(x/b(2)).^b(3));

