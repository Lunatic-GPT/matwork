function [fa,s,pwr,err]=TSE_angles_Busse(Iss,necho,T2,T1,tau,fa_exc,necho_const,no_rampdown)
% [fa,s,pwr,err]=TSE_angles_Busse(Iss,necho,T2,T1,tau,fa_exc,necho_const)
% output:
% fa: angles of the 
% s: the signal intensity time course 
% pwr: pwr level in units of rad^2
% Input:
% Iss: steady state signal
% necho: total echo train length
% tau: half echo spacing
% fa_exc: the first excitatoin pulse angle in degrees
% necho_const: Number of echos that maintain constant signal intensity

if ~exist('no_rampdown')
    no_rampdown = false;
end

pwr=-1;

n=necho*2+2;
f=zeros(length(fa_exc),2*n+1);
z=zeros(length(fa_exc),2*n+1);


%precess
f(:,n+1)=sin(fa_exc*pi/180);  
z(:,n+1)=cos(fa_exc*pi/180);  % ... f(n-2), f(n), f(n+2) ...; ... z(n-1), z(n+1), z(n+3), ... do not contribute to f(n+1);

%but z(n-1), z(n+1), z(n+3) do contribute to the signal of the next TR

s=[];
fa=[];
za=zeros(10,2*n+1);

Issa=zeros(1,necho);
if no_rampdown
    Issa(:)=Iss;
else

for i=1:necho
    if i==1
        Issa(i)=(1+Iss)/2;
    else
        Issa(i)=(Iss+Issa(i-1))/2;
    end
end

end

for irep=1:necho_const
    
z2=zeros(length(fa_exc),2*n+1);
f2=zeros(length(fa_exc),2*n+1);
for i=1:2*n+1
    if i>1
    f2(:,i)=f(:,i-1)*exp(-tau/T2);
    end
    if i~=n+1
       z2(:,i)=z(:,i)*exp(-tau/T1);
    else
       z2(:,n+1)=1-(1-z(:,n+1))*exp(-tau/T1);  
    end
end

z=z2;
f=f2;
% 
% if irep==1
%     a=pi;
% elseif irep<=necho_ramp
   [a,err]=desired_angle_TSE(z(n+2),f(n+2),f(n),Issa(irep)*exp(tau/T2));
  
   if err 
       disp('An error occurred during calculation');
       return;
   end
   
    a=a(2);
% else
%     a=desired_angle_TSE(z(n+2),f(n+2),f(n),Issa(irep));
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
    if i~=n+1
       z2(:,i)=z(:,i)*exp(-tau/T1);
    else
       z2(:,n+1)=1-(1-z(:,n+1))*exp(-tau/T1);  
    end
    
end

z=z2;
f=f2;

s(:,irep)=f(:,n+1);
fa(:,irep)=a*180/pi;
za(irep,:)=z2;
end


slp=(fa(end)-fa(end-2))/2;


fa(end+1:necho)=fa(end)+slp*(1:necho-necho_const);

pwr=sum(fa.^2)+fa_exc^2;  

pwr=pwr*(pi/180)^2;

s=calc_signal_TSE_EPG(fa,fa_exc,tau,T1,T2);


