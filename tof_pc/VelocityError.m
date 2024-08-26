function sig_v=VelocityError(sigma,img,venc)
% sig_v=VelocityError(sigma,img,venc)
% sigma is the background noise (std) at regions with no signal
% img: image intensity
% venc: 
% 


SNR=img/sigma*sqrt(2-pi/2);

sig_ph=1/SNR;

sig_v=sig_ph/pi*venc*sqrt(2); % subtraction of two images introduce a factor of sqrt(2);



