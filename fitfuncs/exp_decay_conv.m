function y = exp_decay_conv(b,input,TR)
% y = exp_decay_conv(b,input,TR)
% b(2): time constant
% b(1): the amplitude of the convoluted function 
% b(3): delay
% the input should have positive values.

x=(0:length(input)-1)*TR;

x=exp(-(x-b(3))/b(2)).*(x>=b(3));

x=reshape(x,size(input));
y=conv(input,x);
y=b(1)*y/max(y);

y=y(1:length(input));



