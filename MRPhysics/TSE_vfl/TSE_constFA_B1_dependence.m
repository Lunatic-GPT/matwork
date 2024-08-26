function TSE_constFA_B1_dependence()

%% Alsop's figure; no T1 and T2 relaxations
T1=1.6;
T2=0.08;
tau=0.02:0.02:0.16;  %half the ESP

%%fixed schedule of flip angles

TE=0.153;

necho_const=61;
necho=121;
aexc=90;

B1scl=linspace(0.5,1.5,11);

fa=180*ones(1,2);
%% single echo sequence; changing tau
for i=1:length(B1scl)
  for j=1:length(tau)
    s(:,j,i)=calc_signal(fa*B1scl(i),90*B1scl(i),tau(j),T1,T2);
  end
  
  
end
figure;plot(tau*2,squeeze(s(1,:,:)));
set(gca,'Yscale','log');

%% multiecho sequence; 

B1scl=linspace(0.5,1,6);
tau=0.03;
fa=180*ones(1,6);

for i=1:length(B1scl)

    s2(:,i)=calc_signal(fa*B1scl(i),90*B1scl(i),tau,T1,T2);
  
  
end

figure;plot(tau*2*(1:length(fa)),s2);

set(gca,'Yscale','log');
legend('1','2','3','4','5','6');

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

