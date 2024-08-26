function y = exp_decay_conv_bl(b,input,TR,range)
% y = exp_decay_conv(b,input,TR)
% b(2): time constant
% b(1): the amplitude of the convoluted function 
% b(3): delay
% the input should have positive values.

x=(0:length(input)-1)*TR;

delay=b(5);
x=exp(-(x-delay)/b(2)).*(x>=delay);
%x=exp(-x/b(2));

x=reshape(x,size(input));
y=conv(input,x);
y=b(1)*y/max(y);

y=y(1:length(input));

y=b(1)*y/max(y)+b(3)+b(4)*(x-length(input)/2)/TR/length(input)*2;

y=y(range);