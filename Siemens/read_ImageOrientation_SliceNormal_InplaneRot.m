function read_ImageOrientation_SliceNormal_InplaneRot(fname)
a=readdPar(fname,'ImageOrientationPatient');
extp(fname);
x=readsPar([fname,'.pro'],'sSliceArray.asSlice[0].sNormal.dSag');
y=readsPar([fname,'.pro'],'sSliceArray.asSlice[0].sNormal.dCor');
z=readsPar([fname,'.pro'],'sSliceArray.asSlice[0].sNormal.dTra');
rot=readsPar([fname,'.pro'],'sSliceArray.asSlice[0].dInPlaneRot');

fprintf('ImageOrientationPatient: %s\n',num2str(a'));
fprintf('normal: %s\n',num2str([x,y,z]));
fprintf('Inplane rot: %f\n',rot);



