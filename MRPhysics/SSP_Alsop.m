function s=SSP_Alsop(a)
%a is in degrees
a=a*pi/180;
sa=sin(a);
ca=cos(a);


phi=linspace(0,2*pi,1000);
phi=phi(1:end-1);

y=1./sqrt(1+(sa/(1-ca))^2/2*(1-cos(phi)));

s=sum(y)*(phi(2)-phi(1))/2/pi;
