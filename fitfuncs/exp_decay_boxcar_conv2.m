function y=exp_decay_boxcar_conv2(b,x)
% y=b1*exp(-x/b2);
% %4.3f*exp(-x/%4.3f)
y=b(1)*(1-exp(-(x-b(3))/b(2))).*double(x>b(3))-b(1)*(1-exp(-(x-b(4))/b(2))).*double(x>b(4));

