function TSE_VarFA_CNR
T1=[1400,4200];
T2=[76,790];

npe=round(404/3)+20;

esp=8.62;
T1sim=4200;
T2sim=790;

TE = linspace(200,800,51);
Iss=linspace(0.1,0.5,41);

a=zeros(2,length(TE),length(Iss));
pwr=zeros(length(TE),length(Iss));
necho_const=zeros(1,length(TE));
 matlabpool(3);

parfor ii=1:length(Iss)
    a2=zeros(2,length(TE));
    pwr2=zeros(1,length(TE));
  for st=1:length(TE)
  disp([ii,st]);
          
   tau=esp/2;  %half the ESP
   necho_const=round(TE(st)/tau/2);
   necho=npe/2+necho_const-1;
   aexc=90;
   fa=TSE_angles_Busse(Iss(ii),necho,T2sim,T1sim,tau,aexc,necho_const);
   if length(fa)~=necho || fa(end)>180
     continue;
   end
   
   for i=1:length(T1)
     aexca=90;
   %  [ss,pwr(i,st,ii)]=calc_signal(fa*aexca/90,aexca,tau,T1(i),T2(i));
     [ss,pwr2(st)]=calc_signal(fa*aexca/90,aexca,tau,T1(i),T2(i));
     a2(i,st)=ss(necho_const);
   end
  end
 a(:,:,ii)=a2;
 pwr(:,ii)=pwr2;
end

save TSE_VarFA_CNR_pwr_PVST1T2

%{


TE = [537,707,609,859,859];
Iss =[0.42,0.40,0.36,0.36,0.3];

necho=[277,156,218,248,248];
%npe=248;

esp=[3.96,8.62,8.34,8.34,8.34];


a=zeros(npe,2,length(TE));

for st=1:length(TE)


for i=1:length(T1)
    
    [a(:,i),necho_const] = TSE_VarFA_Signal(T1(i),T2(i),TE(st),Iss(st),npe(st),esp(st));
     
end

fprintf('contrast = %f\n', diff(a(necho_const,:)));

end
%}

function [s,necho_const,pwr,ss]=TSE_VarFA_Signal(T1,T2,TE,fa,npe,esp)

%T1=3.7;
%T2=0.6;
tau=esp/2;  %half the ESP

necho_const=round(TE/tau/2);
necho=npe/2+necho_const;


%fa=TSE_angles_Busse(Iss,necho,600,3500,tau,aexc,necho_const);
%fa=load('FlipAngle_pss_1024.log');
if length(fa)~=necho || fa(end)>180
    s=0;
    pwr=0;
    return;
end

aexca=90;%linspace(45,135,3);
%aexca=135;
%fa2a=zeros(1,30);
 %fa2=TSE_angles_Busse(0.12,necho,T2,T1,tau,aexc,necho_const);
 %figure;
 %sym={'r','b','k'};
     
  [s,pwr]=calc_signal(fa*aexca/90,aexca,tau,T1,T2);
  
  ss=s;
     s=s(necho_const);
 
 %hold on;plot(s,sym{i});



%s=calc_signal(fa,90,tau,T1,T2);

%s2=calc_signal(fa2,90,tau,T1,T2);



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



function [s,pwr,mz]=calc_signal(fa,fa_exc,tau,T1,T2)
%tau is half of echo spacing


fa=fa*pi/180;
aexc=fa_exc*pi/180;
pwr=sum(fa.^2)+aexc^2;

%rfp=(pi^2*length(fa)+pi^2/4);
%pwr=pwr/rfp;


s=[];
mz=[];

n=length(fa)*2+2;

f=zeros(length(aexc),2*n+1);
z=zeros(length(aexc),2*n+1);


%precess
f(:,n+1)=sin(aexc);
z(:,n+1)=1;%cos(aexc);
for irep=1:length(fa)
    
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
mz(irep)=1+z2(n+1);

end

%pwr=pw90*B90^2*ints90+pw180*B180^2*ints180*(sum(aa(b1ind,:,1).^2,2))/(pi^2);

%fprintf('power = %f \n',pwr);
