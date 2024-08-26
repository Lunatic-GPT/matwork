function y=inv_linear(b,x)
% y=b1*(x+b2);
% %4.3f/(x/%4.3f+1)
y=b(1)./((x/b(2)).^2+1);

