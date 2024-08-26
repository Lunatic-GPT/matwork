function [r,a]=MT_rfs(ecc)
%[r,a]=MT_rfs(ecc)
a=0.64*ecc+18.93;  % Ref: Felleman, el.al, Journal of Physiology, vol 52, pp. 488, 1984.  //owl monkey.
r = sqrt(a/pi);

r = 1.04+0.61*ecc;  % albright 
r = 2.1*(ecc.^0.47); % Maunsell
