function cbv = dR2stoCBV(dR2s,R1,df,Hct)
%dR2stoCBV(dR2s,R1) or dR2stoCBV(dR2s,[],df,Hct)
% calculate the blood volume fraction based on R2star changes due to MION injection and R1
% or df.
% R1 is the 1/T1 of blood plasma containing MION.
% df is the frequency difference of the proton in cross tubes caused by
% the plasma.

if ~exist('Hct','var')
    Hct=0.4;
end
if ~isempty(R1)
    df = 116.9*R1-207.1; 
end

cbv=3/4/pi*dR2s/(df*(1-Hct));







