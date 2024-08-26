function y=parabolic_fit(b,x)
% y=b1*(x-b2)^2+b3
% y = %5.4f*(x-%5.4f)^2+%5.4f
y=b(1)*(x-b(2)).^2+b(3);