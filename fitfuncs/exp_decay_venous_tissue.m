function y=exp_decay_venous_tissue(b,x)
% exp_decay_venous_tissue(b,x);
% x is echo time in ms.
% b1: O2 saturation level
% b2: correlation time
% b3: blood volume fraction
% b3*exp(-R_2v*t)+(1-b3)exp(-R_2t*t).
% 
% Y=%3.2f, \tau=%3.2f ms, \nu=%3.2f

 R2v=24+1125*(1-b(2))*(1-2*b(2)./x.*tanh(x./2/b(2)))./(1-2*b(2)/40*tanh(40/2/b(2)));
 
 R2t = 25;
 
 y=b(3)*exp(-R2v.*x)+(1-b(3))*exp(-R2t.*x);
 
