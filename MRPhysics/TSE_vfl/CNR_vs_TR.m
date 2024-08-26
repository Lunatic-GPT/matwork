TR = linspace(2,10);
T1 = [3.5,1.3];
water=[1,0.7];
T2=[600,70];
TE = 300;

s=zeros(length(TR),length(T1));

for i=1:length(TR)
    for j=1:length(T1)
        
     s(i,j)=(1-exp(-(TR(i)-1)/T1(j)))*water(j)*exp(-TE/T2(j));
     
         

    end
end


c=s(:,1)-s(:,2);

c_norm=c./sqrt(TR');

figure;plot(TR,c_norm,'o-');