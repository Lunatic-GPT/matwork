function TSE_VarFA()

%% Alsop's figure; no T1 and T2 relaxations
T1=3.5;
T2=0.6;
tau=0.00804/2;  %half the ESP
%tau=0;
%%fixed schedule of flip angles

TE=0.6;

necho_const=round(TE/tau/2);

necho=round(1.3/tau/2);
aexc=90;
%fa=TSE_angles_TRAPS_HennigMRM03([4,10,10,20],27,60,60);

%fa=TSE_angles_TRAPS(round(TE/tau),tbFactor);

%fa=load('TSE_vfl/FlipAngles_tb183_stb2.log');
%fa2=load('../../ForSiemens/temp/Flip_array_pss.log');
%fa2=TSE_angles_Busse(0.2,183,T2,T1,tau,90,94);
%fa=load('../../ForSiemens/temp/flipangles.log');
%fa=TSE_angles_exp_flat_exp(0.3,necho,T2,T1,tau,aexc,40*tau,120*tau,necho_const);

fa=TSE_angles_Busse(0.4,necho,0.6,T1,tau,aexc,necho_const);
%fa=load('FlipAngle_pss_1024.log');

aexca=linspace(45,135,3);
%aexca=135;
%fa2a=zeros(1,30);
 %fa2=TSE_angles_Busse(0.12,necho,T2,T1,tau,aexc,necho_const);
 figure;
 sym={'r','b','k'};
for i=1:length(aexca)
     
  s=calc_signal(fa*aexca(i)/90,aexca(i),tau,T1,T2);
 hold on;plot(s,sym{i});
end

for i=1:length(aexca)
     
  s=calc_signal(fa2*aexca(i)/90,aexca(i),tau,T1,T2);
 figure;plot(s);
 set(gca,'FontSize',14);
 xlabel('Echo Index');
 ylabel('Intensity (M_0)');
 xlim([0,length(s)]);
end

%fa=[30*ones(1,9),linspace(30,180,31),linspace(175,30,30),30*ones(1,9)];

%fa=TSE_angles_StabAmp(4,0.12,50);
%fa2=[90*ones(1,20),60*ones(1,30)];
%fa=120*ones(1,35);

%B1scl=ones(1,11)';

%for 

s=calc_signal(fa,90,tau,T1,T2);

s2=calc_signal(fa2,90,tau,T1,T2);

[x,ps]=psf(s);

[x2,ps2]=psf(s2);
% 
% TEeff=frac_transverse(fa,90,tau,T1,T2,TE);
% fprintf('TE = %f; effective TE = %f\n',TE,TEeff);
% TEeff=frac_transverse(fa2,90,tau,T1,T2,TE);
% fprintf('TE = %f; effective TE = %f\n',TE,TEeff);
% 

figure;
subplot(2,1,1);
plot(s(:)');

hold on;
plot(s2(:)','r');
%ylim([0,0.1]);
xlim([0,length(fa)]);

%plot(exp(-tau*2*(1:size(s,2))/T2),'k','LineWidth',2);

set(gca,'FontSize',14);
xlabel('Echo Number');
ylabel('Echo Amplitude');
title('Evolution');
subplot(2,1,2);
plot(x,ps);
hold on;

plot(x2,ps2,'r');
xlim([-3,3]);
set(gca,'FontSize',14);
xlabel('Pixel Number');
ylabel('Signal (arb. units)');
title('Point Spread Function');
% text(10,s(6,end)+0.03,'180^o','FontSize',14);
% text(10,s(3,end)+0.03,'126^o','FontSize',14);
% text(10,s(1,end)+0.03,'90^o','FontSize',14);
legend('Siemens','Busse');
figure;plot(fa);
hold on;
plot(fa2,'r');

set(gca,'FontSize',14);
xlabel('Echo Number');
ylabel('Flip Angle (degree)');
%text(10,s(6,end)+0.03,'180^o','FontSize',14);
%text(10,s(3,end)+0.03,'126^o','FontSize',14);
%text(10,s(1,end)+0.03,'90^o','FontSize',14);
xlim([0,length(fa)]);
ylim([0,180]);
legend('Siemens','Busse');

% 
% figure;
% plot(B1scl_ref*180,s(:,TEecho,2),'ko-');
% set(gca,'FontSize',14);
% xlabel('Peak Flip Angle (degree)');
% ylabel('Peak Amplitude');
% 
% xlim([80,280]);

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

