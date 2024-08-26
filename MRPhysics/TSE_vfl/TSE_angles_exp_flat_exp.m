function [aa,err]=TSE_angles_exp_flat_exp(Iss,necho,T2,T1,tau,aexc,TC1,TC2,necho_const)




%
% adEvolutionParameters[sEvolution.ulName][6]: total time for the second
% exp
% adEvolutionParameters[sEvolution.ulName][5]: time const for the second exp
% adEvolutionParameters[sEvolution.ulName][4]: total time for the first exp
% adEvolutionParameters[sEvolution.ulName][3]: time const for the first exp

%{ 0.5, 100.0, 1000.0, 174.7030, 119.4969, 289.8299, 465.4088, 0.0, 0.0, 0.0 }

%
%fa=TSE_angles_exp_flat_exp(Iss,necho,T2,T1,tau,aexc,40*tau,40*tau,necho_const);

if isempty(TC1)
    TC1=0.1747;
end

if isempty(TC2)
    TC2=0.2898;
end

aexc=aexc*pi/180;

n=necho*2+2;
f=zeros(length(aexc),2*n+1);
z=zeros(length(aexc),2*n+1);


%precess
f(:,n+1)=sin(aexc);
z(:,n+1)=cos(aexc);

s=[];
aa=[];
za=zeros(10,2*n+1);

Issa=zeros(1,necho);

for i=1:necho
    if i<=necho_const
      tmp=sin(aexc)*exp(-tau*2*i/TC1);
      if tmp>Iss
        Issa(i)=tmp;
      else
        Issa(i)=Iss;  
      end
    else
        
      Issa(i)=Iss*exp(-tau*2*(i-necho_const)/TC2);
    end
end

%%

for irep=1:necho
    
z2=zeros(length(aexc),2*n+1);
f2=zeros(length(aexc),2*n+1);
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
if irep==1 || Issa(irep) ~= Issa(irep-1) 
    [a,err]=desired_angle_TSE(z(n+2),f(n+2),f(n),Issa(irep));
    if err
        return;
    end
    a=a(2);
else
     a=-2*atan((Iss-f(n))/z(n+2)/2);
end
     %     
% end
%a=min(a);
z2=zeros(length(aexc),2*n+1);
f2=zeros(length(aexc),2*n+1);
for i=1:2*n+1
   mi=2*n+2-i;
   f2(:,i)=(1+cos(a)).*f(:,i)/2+(1-cos(a)).*f(:,mi)/2+sin(a).*z(:,i);
   z2(:,i)=-f(:,i).*sin(a)/2+f(:,mi).*sin(a)/2+cos(a).*z(:,i);
end

z=z2;
f=f2;

z2=zeros(length(aexc),2*n+1);
f2=zeros(length(aexc),2*n+1);
for i=1:2*n+1
    if i>1
     f2(:,i)=f(:,i-1)*exp(-tau/T2);
    end
     z2(:,i)=z(:,i)*exp(-tau/T1);
end

z=z2;
f=f2;

s(:,irep)=f(:,n+1);
aa(:,irep)=a;
za(irep,:)=z2;
end

aa=aa*180/pi;
err=0;
