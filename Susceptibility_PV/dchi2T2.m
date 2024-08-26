function T2=dchi2T2(dchi,Hct,chi_do,chi_o)

% parameters from Blockley et al., Magnetic Resonance in Medicine 60:1313–1320 (2008)
if ~exist('Hct','var')
Hct=0.42; % male;
end

if ~exist('chi_do','var')
    chi_do=0.27*4*pi*1e-6;  % value from Fan et al., MRM, 72:149-159 (2014)
end

if ~exist('chi_o','var')
    chi_o=-0.03*4*pi*1e-6;
end

Yv=dchi2Yv(dchi,Hct,chi_do,chi_o);

R2=60.45*(1-Yv)/0.11*Hct/0.4;
T2=1/R2;