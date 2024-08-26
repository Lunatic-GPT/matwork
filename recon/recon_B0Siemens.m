function recon_B0Siemens(scand,dte)

a=readMeasDat_NavCorr(scand);

fa=fft2c(a);
write_afni(abs(fa(:,:,:,:,1)),[scand,'_recon']);

b=angle(conj(fa(:,:,:,:,1)).*fa(:,:,:,:,2))/dte/2/pi*1000;  %b in hertz.


write_afni(b,[scand,'_B0']);

