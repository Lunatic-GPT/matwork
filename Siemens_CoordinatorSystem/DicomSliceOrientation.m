function [norm,ori]=DicomSliceOrientation(dname)

if exist(dname,'dir')
    f=name4pat(fullfile(dname,'*.IMA'),1);
end
a=dicominfo(f);
o=a.ImageOrientationPatient;
norm=cross(o(1:3),o(4:6));
ori=Norm2SliceOrient(norm);

fprintf('Norm = %s; %s\n',num2str(norm(:)'),ori);
