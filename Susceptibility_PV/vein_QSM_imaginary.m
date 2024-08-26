function [dchi,flag]=vein_QSM_imaginary(S,TE,B0,T2s,theta_B0)

gamma=42.58e6*2*pi;  %rad/T/s


c = (1-3*cos(theta_B0*pi/180)^2)/3*0.5*gamma*B0*1e-6;

ratio=imag(S(1))*exp(-TE(2)/T2s)/imag(S(2))/exp(-TE(1)/T2s);


func=@(x) sin(x*c*TE(1)/1000)./sin(x*c*TE(2)/1000)-ratio;

[dchi,tmp,flag]=fsolve(func,0.44);





