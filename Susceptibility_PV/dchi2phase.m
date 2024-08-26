function [ph,maxph]=dchi2phase(a,dchi,B0,TE,R,theta,phi,epislon)
% [ph,maxph]=dchi2phase(a,dchi,B0,TE,R,theta,phi)
% B0: unit T
% TE: unit ms
% theta: degree
% phi: in degree; angle between the probe vector and the projection of field on to the perp plane. 
% dchi: ppm; SI units
% R: distance from the center of the vessel
% a: radius of the vessels
% results proportional to R^2 and is a function of a/R.
% Results from: Hsieh; MRI, 33:420-436 (2015)
% output in radian.

if ~exist('phi','var')
    phi=0;
end


if ~exist('epislon','var')
    epislon=1e-6;
end

gamma=42.58e6*2*pi;  %rad/T/s
theta=theta/180*pi;
phi=phi/180*pi;

dchi=dchi*1e-6;
TE = TE/1000;

g=0.5*gamma*dchi*B0*TE;

ph=zeros(size(phi));

ph_in= -g/3*(1-3*cos(theta).*cos(theta));  %consistent with Siemens; same sign as the field.

ph(R<a+epislon)=ph_in;

  ph_out = g*(a./R).^2.*cos(2*phi)*sin(theta)^2; %consistent with Siemens; same sign as the field.

ph(R>a+epislon)=ph_out(R>a+epislon);

maxph=g;














