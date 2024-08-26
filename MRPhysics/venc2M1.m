function M1=venc2M1(venc)
% M1 in G/cm*s^2
% venc in cm/s;
% Note: Siemens M1 results in mT/m*us^2

g=4258*2*pi;  %radian/s/G

M1=pi/g/venc;
 

