function aa=TSE_angles_StabAmp(necho_ramp, Iss,necho)


T2=0.075;
tau=0;
T1=1.4;

B1scl=1;
aexc=pi/2*B1scl;

n=necho*2+2;
f=zeros(length(aexc),2*n+1);
z=zeros(length(aexc),2*n+1);


%precess
f(:,n+1)=sin(aexc);
z(:,n+1)=cos(aexc);


s=[];
aa=[];
za=zeros(10,2*n+1);

Issa=linspace(1,Iss,necho_ramp);

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

 if irep==1
     a=pi;
 elseif irep<=necho_ramp
    a=desired_angle_TSE(z(n+2),f(n+2),f(n),Issa(irep));
    a=a(2);
else 
    a=-2*atan((Iss-f(n))/z(n+2)/2);
end
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
%figure;plot(s);

