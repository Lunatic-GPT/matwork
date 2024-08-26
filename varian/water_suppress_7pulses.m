function water_suppress_7pulses()

flip1=60;
db=4;
flip2=flip1*10^(db/20);
flip=[flip1,flip1,flip1+flip2,flip1,flip1+flip2,flip1,flip1+flip2];
d=[0.15,0.08,0.16,0.08,0.1,0.035,0.078];

s=1;
T1 =2;
for i=1:7
   s=s*cos(flip(i)*pi/180);
   
   s=1-(1-s)*exp(-d(i)/T1);
end
disp(s);
