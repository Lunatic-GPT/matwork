function res=pO2_to_SO2(pO2)
% pO2 in mmHg 
% Severinghaus: 1979 Eq. 1
res=100./(1./(pO2.^3+150*pO2)*23400+1);



