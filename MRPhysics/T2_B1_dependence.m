function T2_B1_dependence()

%% Alsop's figure; no T1 and T2 relaxations
T1=3.6;
T2=0.6;

%T1=1.94;
%T2=0.091;
T1=1.40;
T2=0.067;
%T1=4.2;
%T2=0.79;


tau=0.0127/2;  %half the ESP
tbf=12;
%TE = [912,760,608,456,304,76]/1000;
TE=[183,148,113,78,44,8.7]/1000;
tau=0.0087/2;
tbf=4;

%%fixed schedule of flip angles

B1scl=linspace(0.67,1.33,21);

fa=180*ones(tbf*6,1);

T2m=zeros(1,length(B1scl));

for i=1:length(B1scl)

%figure;
    s(:,i)=calc_signal(fa*B1scl(i),90*B1scl(i),tau,T1,T2);

 y=s(round(TE/tau/2),i);
 x=TE;
 options=optimset('MaxFunEvals',40);
 beta=lsqcurvefit(@exp_decay,[max(y),mean(TE)],x,y',[],[],options);

 T2m(i)=beta(2);
 
 plot(x,y,'o');
 hold on;
 plot(x,exp_decay(beta,x),'r-');
 set(gca,'yscale','log');
 
end

figure;plot(B1scl,T2m);

%save T2_B1_dependence_GM T2m B1scl 
save T2_B1_dependence_WMperp T2m B1scl 


function [s,pwr]=calc_signal(fa,fa_exc,tau,T1,T2)
%tau is half of echo spacing


s=[];

aexc=fa_exc*pi/180;
fa=fa*pi/180;
n=400;

f=zeros(length(aexc),2*n+1);
z=zeros(length(aexc),2*n+1);


%precess
f(:,n+1)=sin(aexc);
z(:,n+1)=cos(aexc);
for irep=1:length(fa)
    
%precess    


z2=zeros(1,2*n+1);
f2=zeros(1,2*n+1);
for i=1:2*n+1
    if i>1
    f2(:,i)=f(:,i-1)*exp(-tau/T2);
    end
    z2(:,i)=z(:,i)*exp(-tau/T1);
end

z=z2;
f=f2;

a=fa(irep);
%a=min(a);
z2=zeros(1,2*n+1);
f2=zeros(1,2*n+1);
for i=1:2*n+1
   mi=2*n+2-i;
   f2(:,i)=(1+cos(a)).*f(:,i)/2+(1-cos(a)).*f(:,mi)/2+sin(a).*z(:,i);
   z2(:,i)=-f(:,i).*sin(a)/2+f(:,mi).*sin(a)/2+cos(a).*z(:,i);
end

z=z2;
f=f2;

z2=zeros(1,2*n+1);
f2=zeros(1,2*n+1);
for i=1:2*n+1
    if i>1
     f2(:,i)=f(:,i-1)*exp(-tau/T2);
    end
     z2(:,i)=z(:,i)*exp(-tau/T1);
end

z=z2;
f=f2;

s(irep)=f(:,n+1);

end

%pwr=pw90*B90^2*ints90+pw180*B180^2*ints180*(sum(aa(b1ind,:,1).^2,2))/(pi^2);
pwr=sum(fa.^2)+aexc^2;

rfp=(pi^2*length(fa)+pi^2/4);
pwr=pwr/rfp;
fprintf('power = %f \n',pwr);



%%

