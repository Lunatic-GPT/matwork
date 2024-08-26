function res=radiationDose(dose,halfLife,bodyWeight,absorptionFraction)

% dose in mCi
% halfLife in min
% bodyWeight in kg

T = halfLife*60/log(2);
E_ev=dose*3.7e7*T*511000*2; % in eV
E_J=E_ev/6.2e18;

res=E_J/bodyWeight*absorptionFraction; %Sv
