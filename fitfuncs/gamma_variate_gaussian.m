function y=gamma_variate_gaussian(b,x)
% y = b1*(x-b2)^b3*exp(-(x-b2)/b4)+b5*(x-b6)^b7*exp(-(x-b6)/b8);
% %3.2f(x-%4.3f)^{%3.2f}exp(-x/%3.2f)+%3.2f(x-%4.3f)^{%3.2f}exp(-x/%3.2f)
y = b(1)*(x-b(2)).^b(3).*exp(-(x-b(2))/b(4)).*(x>b(2))+b(5)*exp(-(x-b(6)).^2/b(7)/b(7)/2);