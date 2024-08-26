function recon2D_siemens(fname)

d=readMeasDat_NavCorr(fname);

fd=fft2c(d);

prefix=strtok(fname,'.');
write_afni(abs(fd),prefix);

