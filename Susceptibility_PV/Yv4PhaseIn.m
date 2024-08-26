function [Yv,dchi]=Yv4PhaseIn(ph,TE,theta_B0,B0)
% TE: in ms
% theta_B0: degree
% ph in radian

gamma=42.58e6*2*pi;
TE = TE/1000;
g=ph*3/(3*cos(theta_B0*pi/180)^2-1);

dchi=g*2/gamma/B0/TE;
Yv=dchi2Yv(dchi);

if nargout==0
    fprintf('Chi = %f ppm\n',dchi*1e6);
    fprintf('Yv = %f\n',Yv);
    
end
    

