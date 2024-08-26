function res=mean_angle(d)
% d in degree;

d=d*pi/180;

res=angle(mean(exp(1i*d)))*180/pi;

