fa=linspace(0,30);
M=zeros(1,length(fa));
for i=1:length(fa)
    M(i)=ss_fit([1,fa(i),2]);
end

figure;plot(fa,M,'r-');

M=zeros(1,length(fa));
for i=1:length(fa)
    M(i)=ss_fit([1,fa(i),1.3]);
end


Ea1=Ernst_Angle(2);
Ea2=Ernst_Angle(1.3);
hold on;plot(fa*Ea2/Ea1,M,'k--');

