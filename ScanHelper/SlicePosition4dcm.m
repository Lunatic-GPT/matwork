function SlicePosition4dcm(fname,name_out)

    a=readdPar(fname,'ImageOrientationPatient');
    fprintf('ImageOrientationPatient in dicom header %s\n',num2str(a(:)'));
  extp(fname);
  x=readsPar([fname,'.pro'],'sSliceArray.asSlice[0].sNormal.dSag');
  y=readsPar([fname,'.pro'],'sSliceArray.asSlice[0].sNormal.dCor');
  z=readsPar([fname,'.pro'],'sSliceArray.asSlice[0].sNormal.dTra');
  rot=readsPar([fname,'.pro'],'sSliceArray.asSlice[0].dInPlaneRot');
  
  

center(1)=readsPar([fname,'.pro'],'asSlice[0].sPosition.dSag');
center(2)=readsPar([fname,'.pro'],'asSlice[0].sPosition.dCor');
center(3)=readsPar([fname,'.pro'],'asSlice[0].sPosition.dTra');

  norm=[x,y,z];
  
  fprintf('Read Normal and Inplane rotation angle from dicom: %s - %f deg\n',num2str(norm(:)'),rot*180/pi);

  save_slicePosition(name_out,center,norm,rot);
  
  