function s=dMsatRec_dR1(TR,T1)
%s=dMss_dR1(fa,T1)
% Calculates dM/dR1/M;

if ~exist('T1','var')
    T1=2;
end

s=TR*exp(-TR/T1)/(1-exp(-TR/T1));





