function anesthetic_calc(r_pan,r_ana,r_glc,vol_hr,mass,vol_total,d_ana,d_pan)
% anesthetic_calc(r_pan,r_ana,r_glc,vol_hr,mass,vol_total,d_ana,d_pan)
% r_pan: pancuronium rate. unit: mg/kg/hr.  (assumed concentration:
% 5mg/ml)
% r_ana: anesthetic rate. unit: mg/kg/hr.
% r_glc: glucose rate. default: 0.05 ml/kg/hr.
% vol_hr: infused volume per hour. unit: ml/hr.
% mass: animal mass. unit: kg.
% vol_total: total volume to make.
% d_ana: anesthetic concentration. unit: mg/ml.  
% d_pan: pancuronium concentration. default 5mg/ml
% Alpha Chloralose: 15 mg/ml, 25 mg/ml.
% Dexdormitor: 0.5 mg/ml.

if ~exist('d_pan','var')
    d_pan = 5;
end
d_glc=0.5;  %Dextrose: 0.5 ml/ml.

if isempty(r_glc)
  r_glc = 0.05;  %glucose infusion rate 0.05 ml/kg/hr.
end
% set up equation x*v = b.
% v is [v_pan,v_ana,v_glu,v_sal]'.
x = zeros(4,4);

vm=vol_hr/mass;
x(1,:) = [d_pan*vm/r_pan-1,-1,-1,-1];
x(2,:) = [-1,d_ana*vm/r_ana-1,-1,-1];
x(3,:) = [-1,-1,d_glc*vm/r_glc-1,-1];
x(4,:) = [1,1,1,1];
b=[0;0;0;vol_total];

v = x\b;

fprintf('Pancuranium (%3.1f mg/ml): %3.2f ml\n',d_pan,v(1));
fprintf('Anesthetic (%3.1f mg/ml): %3.2f ml\n',d_ana,v(2));
fprintf('Dextrose (%3.1f%%): %3.2f ml\n',d_glc*100,v(3));
fprintf('Saline: %3.2f ml\n',v(4));












