function res=gradMomentStep(fov)
% fov in mm
% res in mT/m*us
gamma=42.58*2*pi; % MHz/T
fov=fov/1000; %m

res=2*pi/gamma/fov;
res=res*1000;





  