function res = ifft3c(x)

res=ifft1c(ifft1c(ifft1c(x,1),2),3);

