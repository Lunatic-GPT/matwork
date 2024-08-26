clear all;

% change this to your filename
fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves1');

% read the signal in as a complex FID
signal1 = dicom_get_spectrum(fd);

fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves2');
signal2 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves3');
signal3 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves3');
signal3 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves4');
signal4 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves5');
signal5 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves6');
signal6 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves7');
signal7 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves8');
signal8 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves9');
signal9 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves10');
signal10 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves11');
signal11 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves12');
signal12 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves13');
signal13 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves14');
signal14 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves15');
signal15 = dicom_get_spectrum(fd);
fclose(fd);

fd = dicom_open('/data/1/users/gmr/work/academic/phd/data/phantom_fotis/4aves16');
signal16 = dicom_get_spectrum(fd);
fclose(fd);

disp('Done');