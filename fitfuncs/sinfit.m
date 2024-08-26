function y=sinfit(b,x)
%y=cosfit(b,x)
%y=%3.2f*sin(%3.2f*x+%3.2f)+%3.2f;

y=b(1)*sin(b(2)*x+b(3))+b(4);