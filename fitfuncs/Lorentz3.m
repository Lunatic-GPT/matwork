function y=Lorentz3(b,x)
% y=b1/((x/b2)^2+1)+b3;
% %4.3f/((x/%4.3f)^2+1)

y=b(1)./((x/b(2)).^2+1);
