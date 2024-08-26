function y=R1rho(b,x,delta)
%y=Rex(b,x,delta)
% x has two columns:
% 1st column: Omega offset from water in rad/s.  
% 2nd column: omega1 B1 field in rad/s.
%b1 is k
%b2 is p
%b3 is r2
%b4 is r1
sinth2 = x(:,2).^2./(x(:,1).^2+x(:,2).^2);
rex=b(1)*delta^2*b(2)./((delta-x(:,1)).^2+x(:,2).^2+b(1)^2);
y=b(4)*(1-sinth2)+(rex+b(3)).*sinth2;
