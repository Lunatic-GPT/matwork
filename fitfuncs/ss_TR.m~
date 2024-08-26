function s=ss_TR(b,x)
%ss signal as a function of flip angle and T. 
% T1=%4.3ffa=%4.3fM0=%4.3f
fa=b(2);
T1=b(1);
fa=fa*pi/180;
e1=exp(-x/T1);
s=b(3)*sin(fa)*(1-e1)./(1-cos(fa)*e1);

