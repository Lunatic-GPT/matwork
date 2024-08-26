function [Iss_final,fa]=TSE_vfl_At_Fixed_Power(pwri,necho,T2sim,T1sim,tau,aexc,necho_const)
% [Iss_final,fa]=TSE_vfl_At_Fixed_Power(pwr_wanted,necho,T2sim,T1sim,tau,aexc,necho_const)
% pwr is calculated as the sum square of all flip angles (excite and refocus).  Angles are in units of radians.
% tau is half the echo spacing

pwrc=0;
Issc=0;

Isst=1;
Issb=0;
pwrb=0;
pwrt=100000;
while abs(pwrb-pwrt) >0.01
    
    Iss=(Isst+Issb)/2;
    [fa,err,pwr]=TSE_angles_Busse(Iss,necho,T2sim,T1sim,tau,aexc,necho_const);

 if length(fa)~=necho || fa(end)>180
    Isst=Iss;
    Iss=(Iss+Issb)/2;
    pwrt=pwr;
 else
     pwrc(end+1)=pwr;
     Issc(end+1)=Iss;
     
     Issb=Iss;
     Iss=(Iss+Isst)/2;
     pwrb=pwr;
 end
   
end

Iss_final=zeros(1,length(pwri));
fa=[];
for i=1:length(pwri)

  Iss=interp1(pwrc,Issc,pwri(i));

  if isnan(Issc)
      continue;
  end
  
 [fa(:,i),err,pwrt]=TSE_angles_Busse(Iss,necho,T2sim,T1sim,tau,aexc,necho_const);
  
 while abs(pwri(i)-pwrt)>0.01; 
    pwrc(end+1)=pwrt;
    Issc(end+1)=Iss;
    [Issc,ind]=sort(Issc);
    pwrc=pwrc(ind);
    Iss=interp1(pwrc,Issc,pwri(i));
    
  
    [fa(:,i),err,pwrt]=TSE_angles_Busse(Iss,necho,T2sim,T1sim,tau,aexc,necho_const);
 end
 
 Iss_final(i)=Iss;
  
end
