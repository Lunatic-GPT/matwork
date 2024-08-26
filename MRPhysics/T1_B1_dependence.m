function T1_B1_dependence()

%% Alsop's figure; no T1 and T2 relaxations
T1=3.6;
T2=0.6;

T1=[1.93,1.40,4.2];
tau=0.0127/2;
tbf=8;

%%fixed schedule of flip angles

B1scl=linspace(0.5,1.5,31);

T1meas=zeros(length(T1),length(B1scl));

TR = [0.154, 0.500, 1, 2, 4, 8];
for iT1=1:length(T1)
for i=1:length(B1scl)
%figure;

M=zeros(1,length(TR));

 for iTR=1:length(TR)
     sz=1;
     for j=1:15
        
        sz=sz*cos(pi/2*B1scl(i));
        sz=1+(sz-1)*exp(-tau/T1(iT1));
           
        for k=1:tbf
            sz=sz*cos(pi*B1scl(i));
           sz=1+(sz-1)*exp(-tau*2/T1(iT1));
        end
  
        sz=1+(sz-1)*exp(-(TR(iTR)-tbf*tau*2-tau)/T1(iT1));
     end
     
    M(iTR)=sz;
    
 end
 
 y=M;
 x=TR;
 options=optimset('MaxFunEvals',40);
 beta=lsqcurvefit(@exp_decay2,[-max(y),mean(x),max(y)],x,y,[],[],options);

 T1meas(iT1,i)=beta(2);
 
 plot(x,y,'o');
 hold on;
 plot(x,exp_decay2(beta,x),'r-');
 
 
end

 
end

figure;plot(B1scl,T1meas);

%save T2_B1_dependence_GM T2m B1scl 
save T1_B1_dependence T1meas B1scl T1

%%




function [s,pwr,sz]=calc_signal(fa,fa_exc,tau,T1,T2)
%tau is half of echo spacing


s=[];
sz=[];

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
    z2(:,n+1)=1+(z(:,n+1)-1)*exp(-tau/T1);
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
     
    z2(:,n+1)=1+(z(:,n+1)-1)*exp(-tau/T1);
end

z=z2;
f=f2;

s(irep)=f(:,n+1);
sz(irep)=z(:,n+1);

end

%pwr=pw90*B90^2*ints90+pw180*B180^2*ints180*(sum(aa(b1ind,:,1).^2,2))/(pi^2);
pwr=sum(fa.^2)+aexc^2;

rfp=(pi^2*length(fa)+pi^2/4);
pwr=pwr/rfp;
fprintf('power = %f \n',pwr);



%%

