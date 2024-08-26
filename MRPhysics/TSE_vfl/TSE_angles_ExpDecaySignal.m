function [fa,s,pwr,err]=TSE_angles_ExpDecaySignal(S0,T2Sig,necho,T2,T1,tau,fa_exc)
% [fa,s,pwr,err]=TSE_angles_Busse(Iss,necho,T2,T1,tau,fa_exc,necho_const)
% output:
% fa: angles of the 
% s: the signal intensity time course 
% pwr: pwr level; in units of rad^2;
% Input:
% S0: initial signal
% T2sig:  the time constant for signal decay 
% necho: total echo train length
% tau: half echo spacing
% fa_exc: the first excitatoin pulse angle in degrees

pwr=-1;
fa_exc=fa_exc*pi/180;

n=necho*2+2;
f=zeros(length(fa_exc),2*n+1);
z=zeros(length(fa_exc),2*n+1);


%precess
f(:,n+1)=sin(fa_exc);
z(:,n+1)=cos(fa_exc);


s=[];
fa=[];
za=zeros(10,2*n+1);


S_goal=S0*exp(-(0:necho-1)*2*tau/T2Sig);


for irep=1:necho
    
z2=zeros(length(fa_exc),2*n+1);
f2=zeros(length(fa_exc),2*n+1);
for i=1:2*n+1
    if i>1
    f2(:,i)=f(:,i-1)*exp(-tau/T2);
    end
    z2(:,i)=z(:,i)*exp(-tau/T1);
end

z=z2;
f=f2;
% 
% if irep==1
%     a=pi;
% elseif irep<=necho_ramp
   [a,err]=desired_angle_TSE(z(n+2),f(n+2),f(n),S_goal(irep));
   
   if err
       fprintf('An error occurred during calculation at step %d\n',irep);
       
       figure;subplot(2,1,1);plot(fa);
       title('Flip Angle');
       subplot(2,1,2);plot(s,'o');
       hold on;
       plot(S_goal,'r-');
       title('Signal');
       return;
   end
   
    a=a(2);
% else
%     a=desired_angle_TSE(z(n+2),f(n+2),f(n),S_goal(irep));
%     a=a(2);
%     
   % a=-2*atan((Iss-f(n))/z(n+2)/2);
%     
% end
%a=min(a);
z2=zeros(length(fa_exc),2*n+1);
f2=zeros(length(fa_exc),2*n+1);
for i=1:2*n+1
   mi=2*n+2-i;
   f2(:,i)=(1+cos(a)).*f(:,i)/2+(1-cos(a)).*f(:,mi)/2+sin(a).*z(:,i);
   z2(:,i)=-f(:,i).*sin(a)/2+f(:,mi).*sin(a)/2+cos(a).*z(:,i);
end

z=z2;
f=f2;

z2=zeros(length(fa_exc),2*n+1);
f2=zeros(length(fa_exc),2*n+1);
for i=1:2*n+1
    if i>1
     f2(:,i)=f(:,i-1)*exp(-tau/T2);
    end
     z2(:,i)=z(:,i)*exp(-tau/T1);
end

z=z2;
f=f2;

s(:,irep)=f(:,n+1);
fa(:,irep)=a;
za(irep,:)=z2;
end

fa=fa*180/pi;

pwr=sum(fa.^2)+fa_exc^2;  
pwr=pwr*(pi/180)^2;

figure;subplot(2,1,1);plot(fa);
subplot(2,1,2);plot(s,'o');
hold on;
plot(S_goal,'r-');

%figure;plot(s);

