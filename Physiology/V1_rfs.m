function [f,a]=V1_rfs(ecc)

x = log10(ecc*60)-1.5;
y = 1.1438 + 0.1920*x + 0.0714*x^2 + 0.0619*x^3; %Ref: Dows: Exp Brain Res. 1981.
f = (10^y)/60;
a=f*f;