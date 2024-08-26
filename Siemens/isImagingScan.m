
function res = isImagingScan(evalInfo)

b1=determineBitFields(evalInfo,'MDH_REFPHASESTABSCAN');
b2=determineBitFields(evalInfo,'MDH_PHASESTABSCAN');
b3= determineBitFields(evalInfo,'MDH_PHASCOR');
b4=determineBitFields(evalInfo,'MDH_PATREFSCAN');
b5=determineBitFields(evalInfo,'MDH_PATREFANDIMASCAN');
b6=determineBitFields(evalInfo,'MDH_NOISEADJSCAN');
res=~( b1 | b2 | b3 | (b4 & ~b5)  | b6);