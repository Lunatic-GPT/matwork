function dchi=phIn2dchi(phi,B0,TE,theta)
%  dchi=phIn2dchi(phi,B0,TE,theta)
% phi: phase shift in degree.
% B0: unit T
% TE: unit ms
% theta: degree
% dchi: ppm; SI units


gamma=42.58e6*2*pi;  %rad/T/s
theta=theta/180*pi;
TE = TE/1000;
phi=phi/180*pi;
% %%
% dchi=dchi*1e-6;
% 
% g=0.5*gamma*dchi*B0*TE;
% 
% ph=zeros(size(phi));
% 
% ph_in= -g/3*(1-3*cos(theta).*cos(theta));  %consistent with Siemens; same sign as the field.

%%

g=-phi*3/(1-3*cos(theta).*cos(theta));

dchi=g/(0.5*gamma*B0*TE);

dchi=dchi*1e6; % to ppm














