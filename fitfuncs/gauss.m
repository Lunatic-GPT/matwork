function y=gauss(b,x)
%normalized gaussian distribution
y = 1/sqrt(2*pi)/b(2)*exp(-(x-b(1)).^2/2/b(2)/b(2))*b(3);