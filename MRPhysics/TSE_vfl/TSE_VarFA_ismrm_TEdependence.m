function TSE_VarFA_ismrm_TEdependence()

%% Alsop's figure; no T1 and T2 relaxations

T1a=[1.801,1.152,2.647];
T2a=[0.071,0.064,0.343];

T1=T1a(3);
T2=0.343;
tau=0.00434/2;  %half the ESP
%tau=0;
%%fixed schedule of flip angles



necho_const=135;
necho=270;
aexc=60;
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
s=zeros(necho,3);
 for i=1:3    
  [s(:,i),pwr(i)]=calc_signal(fa*0.5,aexca*0.5,tau,T1a(i),T2a(i));
plot(s(:,i),sym{i}); 
hold on;
 end
 
 figure;
 plot(s(:,3)-s(:,1),'r');
 hold on;
 plot(s(:,3)-s(:,2),'b');
 
 set(gca,'FontSize',12);
 xlabel('Echo Index');
 ylabel('S (M_0)');
 legend('GM','WM','CSF');
 title('(A)');
 xlim([0,length(s1)]);
 
 
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


