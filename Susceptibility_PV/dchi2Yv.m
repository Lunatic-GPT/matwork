function Yv=dchi2Yv(dchi,Hct,chi_do,chi_o)

% SI unit
%value from Magnetic Resonance in Medicine 67:669–678 (2012)
%VC 2011
if ~exist('Hct','var')
Hct=0.42; % male;
end

if ~exist('chi_do','var')
    chi_do=0.27*4*pi*1e-6;  % value from Fan et al., MRM, 72:149-159 (2014)
end

if ~exist('chi_o','var')
    chi_o=-0.03*4*pi*1e-6;
end

Yv=1-(dchi-chi_o*Hct)/chi_do/Hct;