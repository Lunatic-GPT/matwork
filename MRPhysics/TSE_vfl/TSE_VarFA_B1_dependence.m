function TSE_VarFA_B1_dependence()

%% Alsop's figure; no T1 and T2 relaxations
T1=1.6;
T2=0.1;
epi_factor =3;

tau=0.00434/2*epi_factor;  %half the ESP

%%fixed schedule of flip angles

TE=0.153;

necho_const=round(61/epi_factor);
necho=round(121/epi_factor);

aexc=90;
%fa=TSE_angles_TRAPS_HennigMRM03([4,10,10,20],27,60,60);

%fa=TSE_angles_TRAPS(round(TE/tau),tbFactor);

%fa=load('TSE_vfl/FlipAngles_tb183_stb2.log');
%fa=load('FlipAngles.log');
%fa=TSE_angles_Busse(0.14,256,T2,T1,tau,90,129);
fa=TSE_angles_exp_flat_exp(0.165,necho,T2,T1,tau,aexc,9*tau,40*tau,necho_const);

fa2=TSE_angles_Busse(0.12,necho,T2,T1,tau,aexc,necho_const);

%fa=[30*ones(1,9),linspace(30,180,31),linspace(175,30,30),30*ones(1,9)];

%fa=TSE_angles_StabAmp(3,0.6,40);
%fa=[90*ones(1,20),60*ones(1,30)];
%fa=120*ones(1,35);

%B1scl=ones(1,11)';

%for 
B1scl=linspace(0.5,1.5,11);


for i=1:length(B1scl)
 s(:,i)=calc_signal(fa*B1scl(i),90*B1scl(i),tau,T1,T2);
 s2(:,i)=calc_signal(fa2*B1scl(i),90*B1scl(i),tau,T1,T2);
 
 
 TEeff(i)=frac_transverse(fa*B1scl(i),90*B1scl(i),tau,T1,T2,TE);

 TEeff2(i)=frac_transverse(fa2*B1scl(i),90*B1scl(i),tau,T1,T2,TE);

end


% 
% TEeff=frac_transverse(fa,90,tau,T1,T2,TE);
% fprintf('TE = %f; effective TE = %f\n',TE,TEeff);
% TEeff=frac_transverse(fa2,90,tau,T1,T2,TE);
% fprintf('TE = %f; effective TE = %f\n',TE,TEeff);

fprintf('TEeff: Siemens Max/Min=%f\n',max(TEeff)/min(TEeff));

fprintf('TEeff: Busse Max/Min=%f\n',max(TEeff2)/min(TEeff2));


fprintf('Signal: Siemens Max/Min=%f\n',max(s(necho_const-1,:))/min(s(necho_const-1,:)));

fprintf('Signal: Busse Max/Min=%f\n',max(s2(necho_const-1,:))/min(s2(necho_const-1,:)));

figure;

plot(B1scl*180,s(necho_const-1,:));
hold on;
plot(B1scl*180,s2(necho_const-1,:),'r');

set(gca,'FontSize',14);
xlabel('Ref. Angle (degrees)');
ylabel('Signal (arb. units)');
%title('Point Spread Function');
legend('Siemens','Busse');


figure;

plot(B1scl*180,TEeff);
hold on;
plot(B1scl*180,TEeff2,'r');

set(gca,'FontSize',14);
xlabel('Ref. Angle (degrees)');
ylabel('Effective TE');
%title('Point Spread Function');
legend('Siemens','Busse');


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

