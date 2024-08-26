function s=dMss_dR1(fa,TR,T1)
%s=dMss_dR1(fa,TR,T1)
% fa in degrees. Calculates dMss/dR1/Mss;

if ~exist('T1','var')
    T1=2;
end


fa2=fa*pi/180;
e1=exp(-TR./T1);
ca=cos(fa2);
sa=sin(fa2);

s=sa*e1*TR/(1-ca*e1)-sa*(1-e1)/(1-ca*e1)^2*ca*e1*TR;

Mss=ss_GEEPI(fa,TR,T1,0);

s=s/Mss;


