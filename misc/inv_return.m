m=1;

r=0.1/12;

N=17;
total=zeros(1,N);

for i=1:N
total(i)=sum(m*(1+r).^(1:i));
end

figure;plot(total./(1:N));

apr=(1+r)^12-1;
disp(apr);