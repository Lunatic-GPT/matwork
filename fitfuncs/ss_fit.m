function s=ss_fit(b,x)
%s=ss_GEEPI(b,x)

if ~exist('x','var')
    x=1;
end
T1=b(3);
TR=0.01;

fa=x*b(2);

fa=fa*pi/180;
e1=exp(-TR/T1);

fa_ernst=Ernst_Angle(TR,T1)*pi/180;
snorm=sin(fa_ernst)*(1-e1)./(1-cos(fa_ernst)*e1);
s=b(1)*sin(fa)*(1-e1)./(1-cos(fa)*e1)/snorm;

