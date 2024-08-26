function [SNR_mag,sig_ph]=imageSNR4backgroundNoise(sigma,img)
% sig_ph is in units of radians
% sigma: std dev of signal intensity in a region devoid of any signal
% img: image intensity 

SNR_mag=img/sigma*sqrt(2-pi/2);  % sqrt(2-pi/2)=0.66;

sig_ph=1./SNR_mag;

