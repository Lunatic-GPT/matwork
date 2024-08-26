function y=gamma_variate(b,x)
% y = b1*(x-b2)^b3*exp(-(x-b2)/b4);
% %3.2f(x-%4.3f)^{%3.2f}exp(-x/%3.2f)

y = b(1)*(x-b(2)).^b(3).*exp(-(x-b(2))/b(4)).*(x>b(2));


% peak at x = b(2)+b(3)*b(4);
% y at peak = 
