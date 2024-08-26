function m=Mss(off,b1,R1,R2)
%m=Mss(off,b1[,R1,R2])
if ~exist('R1','var')
    R1=0.5;
end

if ~exist('R2','var')
    R2=25;
end

A=R2^2/(b1^2+off^2);
th=atan2(b1,off);
m(1)=R1*sin(th)*cos(th)/(R2*sin(th)^2+R1*(cos(th)^2+A));
m(2)=R1*sin(th)*sqrt(A)/(R2*sin(th)^2+R1*(cos(th)^2+A));
m(3)=R1*(cos(th)^2+A)/(R2*sin(th)^2+R1*(cos(th)^2+A));
R1rho=R1*cos(th)^2+R2*sin(th)^2;
m(4)=R1*cos(th)/R1rho*cos(th);  %an approximate solution for m(3)
