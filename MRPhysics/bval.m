function [b,M1]=bval(sep,duration,ampl,trise,moment)
% b=bval(sep,duration,ampl,trise[,moment])
% time in ms
% ampl in G/cm;
% moment: G/cm.ms
% duration is the plateua + one edge
% separation is from center to center
% if duration is empty, then it can be calculated from the moment 
% the final result in units of s/mm2
% M1 is in G/cm*s^2
if ~exist('trise','var')
    trise=0.13*ampl/40;
end

if isempty(duration)
    duration=moment/ampl;
end

trise=trise/1000;
duration=duration/1000;
sep=sep/1000;
g=4258*2*pi;  %1/s/G
b=g^2*ampl^2*duration^2*(sep-duration/3);

b=b+g^2*ampl^2*(trise^3/30-duration*trise^2/6);  %s/cm2

b=b/100;  %s/mm2

M1=duration*sep*ampl;
