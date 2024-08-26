function res=tof_fit(b,x)
% res=tof_fit(b,x)
% 
% b1=%3.2f,b2=%3.2f,b3=%3.2f
res=2*b(1)./(1+exp(-(x/b(2)).^b(3)))+b(4);


