dname='FATNAVIMAGES_20_SPLIT_1_0101';

dir_str=dir(fullfile(dname,'*.IMA'));

a=dicomread(fullfile(dname,dir_str(1).name));

in=dicominfo(fullfile(dname,dir_str(1).name));
in2=dicominfo(fullfile(dname,dir_str(2).name));
 
nslc=length(dir_str); 
x=in.ImageOrientationPatient(1:3);
y=in.ImageOrientationPatient(4:6);
z=cross(x,y);

thk=sos(in2.ImagePositionPatient-in.ImagePositionPatient,1);
in.ImagePositionPatient+in.PixelSpacing(1)*x*floor(size(a,2)/2+0.5)+in.PixelSpacing(2)*y*floor(size(a,1)/2+0.5)-z*(nslc-1)/2*thk



