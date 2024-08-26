dname='16_113037';

dname='22_113037';

disp('TE');
readdPar(dname,'EchoTime')

disp('TR');
readdPar(dname,'RepetitionTime')

disp('dimension');
dcmdim(dname)

disp('FA');
readdPar(dname,'FlipAngle')

disp('TI');
%readdPar(dname,'InversionTime')

disp('Bandwidth');
readdPar(dname,'PixelBandwidth')
extp(dname);

disp('lScanTimeSec');
readsPar([dname,'.pro'],'lScanTimeSec')

disp('thickness');
readsPar([dname,'.pro'],'dThickness')

disp('phaseFOV');
readsPar([dname,'.pro'],'dPhaseFOV')

disp('ReadoutFOV')
readsPar([dname,'.pro'],'dReadoutFOV')

disp('Grappa PE')
readsPar([dname,'.pro'],'lAccelFactPE')

disp('Grappa 3D')
readsPar([dname,'.pro'],'lAccelFact3D')

