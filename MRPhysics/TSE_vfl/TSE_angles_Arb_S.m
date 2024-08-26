function [aa,err,pwr]=TSE_angles_Arb_S(Issa,T2,T1,tau,aexc)

% Issa is the array of signals
% aexc: the excitation RF pulse in degrees

pwr=-1;
aexc=aexc*pi/180;

necho=length(Issa);

n=necho*2+2;
f=zeros(1,2*n+1);
z=zeros(1,2*n+1);

%precess
f(n+1)=sin(aexc);
z(n+1)=cos(aexc);

s=[];
aa=[];
za=zeros(10,2*n+1);

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

    if irep>=25
        disp('');
    end
    [a,err]=desired_angle_TSE(z(n+2),f(n+2),f(n),Issa(irep));
    
    if err
        aa=aa*180/pi;
        return;
    end
    a=a(2);
    
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
    
    s(irep)=f(:,n+1);
    aa(irep)=a;
    za(irep,:)=z2;
end

aa=aa*180/pi;
