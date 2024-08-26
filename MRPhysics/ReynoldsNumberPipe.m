function r = ReynoldsNumberPipe( v,d,rho,visc )
% r = ReynoldsNumberPipe( v,d,rho,visc )
% v: velocity in m/s
% d: dimaeter in m
% rho: density in kg/m3; default 1000 (for water)
% visc: dynamic viscosity in Ns/m2 or equivalently kg/s.m; default: 8.9e-4

if ~exist('visc','var')
    visc=8.9e-4;
end

if ~exist('rho','var')
    rho=1000;
end

r=rho*v*d/visc;

end

