function y=biexp_decay(b,x)
% y=b1*exp(-x/b2);
% %4.3f*exp(-x*%4.3f)+%4.3f*exp(-3.0*x)
y=b(1)*exp(-b(2)*x*0.001)+b(3)*exp(-x*3*0.001);