function plotc(x,y,varargin)

plot(real(x),real(y),varargin{:});
hold on;
plot(imag(x),imag(y),varargin{:});
