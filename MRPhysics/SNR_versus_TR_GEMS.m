T1=2;
TR=linspace(0,3);

for i=1:length(TR)
    
[d,s(i)]=Ernst_Angle(TR(i),T1);

end


figure;plot(TR,s);
hold on;plot(TR,s./sqrt(TR),'r');