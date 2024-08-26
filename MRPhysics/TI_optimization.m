TR = 3.5;
etl=0.962;
TI=linspace(0,2);
T1=1.4;

s0=-exp(-(TR-etl-TI)/T1);

s=1+(s0-1).*exp(-TI/T1);
figure;plot(TI,s);

T1=3.5;

s0=-exp(-(TR-etl-TI)/T1);

s=1+(s0-1).*exp(-TI/T1);
hold on;
plot(TI,s,'r');