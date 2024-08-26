function y=snr_te(b,x)
% y = x2./sqrt(b1^2+(b2*x1^2+b3).*x2.^2);
% x1 is te column vector, x2 is signal intensity column vector.
% b1=sigma0. b2=c1^2*deltaR2^2. b3=c2^2. in Kruger MRM01.

y = x(:,2)./sqrt(b(1)^2+(b(2)*x(:,1).^2+b(3)).*x(:,2).^2);

%y=sin(b(1)*x)+p;