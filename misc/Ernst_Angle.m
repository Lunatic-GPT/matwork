function [d,s]=Ernst_Angle(TR,T1)
% [d,s]=Ernst_Angle(TR[,T1])
% s is the signal at steady state.  s=1 for a pi/2 pulse at infinite TR.  
if ~exist('T1','var')
    T1=2;
end

f=exp(-TR/T1);
d = acos(f);
s=sin(d).*(1-f)./(1-cos(d).*f);

d = d*180/pi;

