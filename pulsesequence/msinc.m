function y=msinc(N,cycle)
% y=msinc(N,cycle)
% cycles on one side (positive or negative).
x = -(N-1)/2:(N-1)/2;
x = x*4*cycle/N;
y = sinc(x);