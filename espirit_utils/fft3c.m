function res = fft3c(x)

res=fft1c(fft1c(fft1c(x,1),2),3);

