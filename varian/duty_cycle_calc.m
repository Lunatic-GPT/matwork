

l=2.56;  %FOV
np=80;
npe=57;
sw=300000;
FOV=2.56;
G=sw/FOV/4258;
TR = 1;
ns=12;

duty_cycle=np/sw*npe/TR*(G/73.4)^2*ns;
disp(duty_cycle);

