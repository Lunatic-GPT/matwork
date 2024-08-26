function [Yv,dchi]=Yv4PhaseOut(ph,TE,theta_B0,phi_B0,B0,r2a)
% [Yv,dchi]=Yv4PhaseOut(ph,TE,theta_B0,phi_B0,B0,r2a)
% TE: in ms
% theta_B0,phi_B0: degree; phi_B0 is the angle of r vec wrt the projection of B0 onto the plane perp to the vessel.   
% ph: radian
% r2a: ratio between distance to center and radius

gamma=42.58e6*2*pi;
TE = TE/1000;

 c=1/r2a^2*cos(2*phi_B0*pi/180)*sin(theta_B0*pi/180)^2;
 
   
    
   g=ph/c;
   
dchi=g*2/gamma/B0/TE;
Yv=dchi2Yv(dchi);



