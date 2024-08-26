function hypercapnia_calc(co2,o2,flow)
% hypercapnia_calc(co2,o2,flow)
x=[0.21-o2,1-o2,0.21-o2; 0.1-co2,-co2,-co2;1,1,1];

b=[0,0,flow]';

c=x\b;

fprintf('10%% CO2 = %f\n O2 = %f\n air = %f\n',c(1),c(2),c(3));


