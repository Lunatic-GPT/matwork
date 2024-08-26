function res=SO2_to_pO2(S)
% pO2 in mmHg 
res=-0.385*log(1/S-1)+3.32-1/(72*S)-0.17*S.^6;
res=exp(res);




