function y=exp_decay_boxcar_conv(b,x)
% y=b1*exp(-x/b2);
% assuming stimulus start from 0
% %4.3f*exp(-x/%4.3f)
y=b(1)*(1-exp(-x/b(2)))-b(1)*(1-exp(-(x-b(3))/b(2)))*int16(x>b(3));

