function res=gradFirstMoment(amp,t,tref)
% amp n*1; 
% t: n*2; the start and end time of the gradients
% a linear ramp is assumed between amplitude transitions
% 
if ~exist('tref','var')
    tref=0;
end

t=t-tref;

res=0;
for i=1:length(amp)
 res=res+amp(i)/2*(t(i,1)^2-t(i,2)^2);
end

for i=1:length(amp)-1  % ramps
    
   t1=t(i,2);
   t2=t(i+1,1);
   k=(amp(i+1)-amp(i))/(t2-t1);
   
   b=amp(i)-k*t1;
   
   res=res+k*(t2^3-t1^3)/3+b*(t2^2-t1^2)/2;
    
    
end
