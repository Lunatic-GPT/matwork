function Signal_vs_FirstFA()

%% the following simulation suggests that only the transverse component after fa_exc matters to the signal intensity.




T1=2.2;
T2=0.5;

tau=0.00434;  %half the ESP

%%fixed schedule of flip angles

TE=0.153;

necho_const=60;
necho=170;

aexc=90;

fa=TSE_angles_Busse(0.12,necho,T2,T1,tau,aexc,necho_const);

B1scl=linspace(0.5,1,20);

s=zeros(length(fa),length(B1scl));

for i=1:length(B1scl)
 s(:,i)=calc_signal_TSE_EPG(fa,90*B1scl(i),tau,T1,T2); 
end


figure;
plot(s);
%%
s2=zeros(length(fa),length(B1scl));

for i=1:length(B1scl)
 s2(:,i)=calc_signal_TSE_EPG(fa,[B1scl(i),0],tau,T1,T2); 
end


figure;
plot(s2);
%%
s3=zeros(length(fa),length(B1scl));

for i=1:length(B1scl)
 s3(:,i)=calc_signal_TSE_EPG(fa,B1scl(i)*[sin(pi/4),cos(pi/4)],tau,T1,T2); 
end

%%


 stemp=calc_signal_TSE_EPG(fa*0.7,[0,1],tau,T1,T2); 

 figure;plot(stemp);

%%
figure;plot(sin(pi/2*B1scl),s(50,:),'o');

hold on;plot(B1scl,s2(50,:),'r');

plot(sin(pi/4)*B1scl,s3(50,:),'ko');





