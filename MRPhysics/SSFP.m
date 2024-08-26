function [s1,s2] = SSFP(theta,T1,T2,TR)
%[s1,s2] = SSFP(theta,T1,T2[,TR])
% s1: FID
% s2:Echo
% default TR = 14 (ms).
if ~exist('TR','var')
    TR = 14;
end

theta = theta.*180/pi;
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);
p=1-E1.*cos(theta)-E2.^2.*(E1-cos(theta));
q=E2.*(1-E1).*(1+cos(theta));

s1=tan(theta/2).*(1-(E1-cos(theta)).*(1-E2.^2)./sqrt(p.^2-q.^2));
s2=tan(theta/2).*(1-(1-E1.*cos(theta)).*(1-E2.^2)./sqrt(p.^2-q.^2));


