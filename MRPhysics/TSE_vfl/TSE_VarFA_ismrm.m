function TSE_VarFA_ismrm()

%% Alsop's figure; no T1 and T2 relaxations

T1a=[1.801,1.152,2.647];
T2a=[0.071,0.064,0.343];

T1=T1a(3);
T2=T2a(3);
tau=0.00434/2;  %half the ESP
%tau=0;
%%fixed schedule of flip angles

TE=0.153;

necho_const=135;
necho=270;
aexc=90;
aexca=90;
%fa=TSE_angles_TRAPS_HennigMRM03([4,10,10,20],27,60,60);

%fa=TSE_angles_TRAPS(round(TE/tau),tbFactor);

%fa=load('TSE_vfl/FlipAngles_tb183_stb2.log');
%fa2=load('../../ForSiemens/temp/Flip_array_pss.log');
%fa2=TSE_angles_Busse(0.2,183,T2,T1,tau,90,94);
%fa=load('../../ForSiemens/temp/flipangles.log');
fa=TSE_angles_exp_flat_exp(0.2,necho,T2,T1,tau,aexc,40*tau,120*tau,necho_const);

%fa2=TSE_angles_Busse(0.12,necho,T2,T1,tau,aexc,necho_const);

%fa2a=zeros(1,30);
 %fa2=TSE_angles_Busse(0.12,necho,T2,T1,tau,aexc,necho_const);
 
figure;
sym={'r','b','k'};
subplot(1,3,1);
 for i=1:3    
  [s1,pwr(i)]=calc_signal(fa*aexca/90,aexca,tau,T1a(i),T2a(i));
  hold on;
  plot(s1,sym{i});
 end
 
 set(gca,'FontSize',12);
 xlabel('Echo Index');
 ylabel('S (M_0)');
 legend('GM','WM','CSF');
 title('(A)');
 xlim([0,length(s1)]);
%%
 subplot(1,3,2);

 
T1=T1a(3);
T2=[0.15,0.155,0.165,0.175,0.185,0.2,0.225,0.25,0.3,0.4,0.6];
tau=0.00434/2;  %half the ESP
%tau=0;
%%fixed schedule of flip angles

necho_const=135;
necho=270;
aexc=90;
aexca=90;
%fa=TSE_angles_TRAPS_HennigMRM03([4,10,10,20],27,60,60);

%fa=TSE_angles_TRAPS(round(TE/tau),tbFactor);

%fa=load('TSE_vfl/FlipAngles_tb183_stb2.log');
%fa2=load('../../ForSiemens/temp/Flip_array_pss.log');
%fa2=TSE_angles_Busse(0.2,183,T2,T1,tau,90,94);
%fa=load('../../ForSiemens/temp/flipangles.log');
s=zeros(necho,3,length(T2));
pwr=zeros(3,length(T2));
for j=1:length(T2)
fa=TSE_angles_exp_flat_exp(0.2,necho,T2(j),T1,tau,aexc,40*tau,120*tau,necho_const);

%fa2=TSE_angles_Busse(0.12,necho,T2,T1,tau,aexc,necho_const);

%fa2a=zeros(1,30);
 %fa2=TSE_angles_Busse(0.12,necho,T2,T1,tau,aexc,necho_const);
 
 for i=1:3    
  [s(:,i,j),pwr(i,j)]=calc_signal(fa*aexca/90,aexca,tau,T1a(i),T2a(i));
 end
end


 plot(T2,squeeze(s(necho/2,3,:)-s(necho/2,1,:)),sym{1});
 hold on;
 plot(T2,squeeze(s(necho/2,3,:)-s(necho/2,2,:)),sym{2});
 
 set(gca,'FontSize',12);
 xlabel('T_2 (ms)');
 ylabel('\Delta S_{TE} (M_0)');
 legend('GM vs CSF','WM vs CSF');
 xlim([0.15,0.6]);
 ylim([0.1,0.3]);
 title('(B)');
 
 
%  %%
%  subplot(1,4,3);
%  
%  T1=T1a(3);
% T2=T2a(3);
% tau=0.00434/2;  %half the ESP
% %tau=0;
% %%fixed schedule of flip angles
% Iss=[0.05,0.12,0.19,0.26,0.295,0.33,0.355,0.365,0.38,0.39,0.4];
% necho_const=135;
% necho=270;
% aexc=90;
% aexca=90;
% %fa=TSE_angles_TRAPS_HennigMRM03([4,10,10,20],27,60,60);
% 
% %fa=TSE_angles_TRAPS(round(TE/tau),tbFactor);
% 
% %fa=load('TSE_vfl/FlipAngles_tb183_stb2.log');
% %fa2=load('../../ForSiemens/temp/Flip_array_pss.log');
% %fa2=TSE_angles_Busse(0.2,183,T2,T1,tau,90,94);
% %fa=load('../../ForSiemens/temp/flipangles.log');
% s2=zeros(necho,3,length(Iss));
% pwr2=zeros(3,length(Iss));
% for j=1:length(Iss)
% fa=TSE_angles_exp_flat_exp(Iss(j),necho,T2,T1,tau,aexc,40*tau,120*tau,necho_const);
% 
% %fa2=TSE_angles_Busse(0.12,necho,T2,T1,tau,aexc,necho_const);
% 
% %fa2a=zeros(1,30);
%  %fa2=TSE_angles_Busse(0.12,necho,T2,T1,tau,aexc,necho_const);
%  
%  for i=1:3    
%   [s2(:,i,j),pwr2(i,j)]=calc_signal(fa*aexca/90,aexca,tau,T1a(i),T2a(i));
%  end
% end
% 
%  plot(Iss,squeeze(s2(necho/2,3,:)-s2(necho/2,1,:)),sym{1});
%  hold on;
%  plot(Iss,squeeze(s2(necho/2,3,:)-s2(necho/2,2,:)),sym{2});
%   set(gca,'FontSize',12);
%  xlabel('S at Plateau (M_0)');
%  ylabel('\Delta S_{TE} (M_0)');
%  legend('GM vs CSF','WM vs CSF');
%  xlim([0.05,0.4]);
%  ylim([0,0.38]);
%  title('(C)');
 
 
 subplot(1,3,3);
 
 plot(squeeze(s(necho/2,3,:)-s(necho/2,2,:)),pwr(1,:),'k-');
 hold on;
 %plot(squeeze(s2(necho/2,3,:)-s2(necho/2,2,:)),pwr2(1,:),'b-');
 set(gca,'FontSize',12);
 xlabel('\Delta S_{TE} (M_0)');
 ylabel('Relative RF Power');
 %legend('Fixed S_{Plateau}','Fixed T_2');
 xlim([0.05,0.35]);
 ylim([0,0.22]);
 title('(C)');
 
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

