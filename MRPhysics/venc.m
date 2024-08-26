function venc=venc(sep,duration,ampl)
% duration in units of ms;
% sep in units of ms;
% ampl in units of G/cm;
% result in cm/s;

 
duration=duration/1000;
sep=sep/1000;
g=4258*2*pi;  %radian/s/G

venc_1 = g*duration*ampl*sep/pi;

venc=1/venc_1;


