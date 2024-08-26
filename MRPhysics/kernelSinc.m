function res=kernelSinc(kmax,x)

res=sinc(kmax*x/pi)*kmax/pi;  % integral over x is 1