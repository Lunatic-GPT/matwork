clear all;


% Phillips' Demo
fd = dicom_open('/home/gmr/work/academic/phd/data/nottingham/DICOM/XX_0106');
y_ph = dicom_get_spectrum_phillips(fd);
fclose(fd);

% Siemens' Demo
fd = dicom_open('/home/gmr/work/academic/phd/data/phantom_single/03311921/16167579');
y_s = dicom_get_spectrum_siemens(fd);
fclose(fd);


disp('done');