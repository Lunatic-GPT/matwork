function s=calc_signal_TSE_EPG(fa,fa_exc,tau,T1,T2)
% calc_signal_TSE_EPG(fa,fa_exc,tau,T1,T2)
% fa_exc: can be one component or two components.  When two components, it gives the transverse and longitudinal components of the initial magnetization. 
% tau is half of echo spacing


s=[];

fa=fa*pi/180;
n=400;


f=zeros(1,2*n+1);
z=zeros(1,2*n+1);

%precess

if length(fa_exc)==1
    
    aexc=fa_exc*pi/180;

    f(:,n+1)=sin(aexc);
    z(:,n+1)=cos(aexc);
else
    f(:,n+1)=fa_exc(1);
    z(:,n+1)=fa_exc(2);
end



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
