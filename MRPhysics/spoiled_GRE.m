function s=spoiled_GRE(theta,TR,T1) 
%s=spoiled_GRE(theta,TR,T1) 
% steady state signal for spoiled GRE
theta = theta/180*pi;
r = TR./T1;
s =sin(theta)*(1-exp(-r))./(1-cos(theta)*exp(-r));


