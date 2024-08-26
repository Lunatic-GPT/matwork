function m=waterDSfit(b,df)

%b1 in Gauss
%df in Hz; can be an array

T1=2;
%T2=0.2;
b1=3e-2;
t=0.5;
m=waterDS(T1,b(2),df-b(1),b1,t);
m=m(:,3)';

m=m*b(3);
      




