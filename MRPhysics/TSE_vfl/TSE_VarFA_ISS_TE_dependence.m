function TSE_VarFA_ISS_TE_dependence()

%% Alsop's figure; no T1 and T2 relaxations

T1a=[1.801,1.152,2.647];
T2a=[0.071,0.064,0.343];

T1=T1a(3);
T2=T2a(3);
tau=0.00396/2;  %half the ESP

TE=[0.3,0.45,0.55,0.695];
Iss=0.05:0.002:0.8;

 tbf=169;
 stb=2;
 
sTE=zeros(3,length(Iss),length(TE));
pwr=zeros(3,length(Iss),length(TE));
s1=cell(3,length(Iss),length(TE));
for k=1:length(TE)
    
for j=1:length(Iss)
    
  necho=round(TE(k)/tau/2)+tbf*stb/2-stb;
  necho_const=round(TE(k)/tau/2);
  
  %[fa,err]=TSE_angles_Busse(Iss(j),necho,T2,T1,tau,90,necho_const);   
[fa,err]=TSE_angles_exp_flat_exp(Iss(j),necho,T2,T1,tau,90,[],[],necho_const);
  if err
      continue;
  end
  
 for i=1:3     
  [s1{i,j,k},pwr(i,j,k)]=calc_signal(fa,90,tau,T1a(i),T2a(i));
 
   sTE(i,j,k)=s1{i,j,k}(necho_const);
   disp([i,j,k]);
 end
 
end

end

save TSE_VarFA_ISS_TE_dependence_siemens

%%
clear all
load TSE_VarFA_ISS_TE_dependence_siemens
figure;
clr={'r','g','b','k'};
for k=1:length(TE)
    ind=(pwr(1,:,k)~=0);
plot(pwr(1,ind,k),sTE(3,ind,k)-sTE(1,ind,k),clr{k});
hold on;
end
 

function TEeff=frac_transverse(fa,fa_exc,tau,T1,T2,TE)

ind=round(TE/tau);
s=calc_signal(fa(1:ind),fa_exc,tau,T1,T2);

s0=calc_signal(fa(1:ind),fa_exc,0,T1,T2);

TEeff=log(s0(end)/s(end))*T2;

function [x,fs2]=psf(s)

n=length(s);
s2=cat(1,s(:),zeros(9*n,1));
fs2=fft1c(s2(:),1);
fs2=abs(fs2);
x=linspace(-(length(s)/2)+0.5,length(s)/2-0.5,length(s)*10);



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

