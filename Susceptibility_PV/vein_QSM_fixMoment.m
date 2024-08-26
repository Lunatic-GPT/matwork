function [res,flag]=vein_QSM_fixMoment(S,rho0,TE,R,B0,theta_B0,dchi_a2)
% S: signal 
% rho0: signal per unit area; unit same as [R]^2
% TE: in ms
% theta_B0 in degrees
% S: signal intensity*[R]^2
% R: in units of mm
% the small angle approach in Hsieh_MRI2015b
% due to digitization error of the total signal, this method give poor
% results when theta_B0>
% dchi_a2: dchi_a2; dchi in units of ppm, a: same units as R


gamma=42.58e6*2*pi; 


fitfunc=@(x) main_func(S,gamma,x(1)*1e-6,B0,TE,theta_B0,rho0,sqrt(dchi_a2/x(1)),R);
 
for i=1:length(TE)
  p_in(i) = dchi2phase(1,0.5,B0,TE(i),0,theta_B0,0);

  if imag(S(i))<=0
      
     lu(i)=-pi/p_in(i)*0.5;
     ll(i)=0;
     
  else
      
     lu(i)=-2*pi/p_in(i)*0.5;
     ll(i)=-pi/p_in(i)*0.5;
      
  end
  
  
end
 
[res,resnorm,tmp2,flag] =      lsqnonlin(fitfunc, (max(ll)+min(lu))/2,max(ll),min(lu),optimoptions('lsqnonlin','TolFun',1e-12,'MaxFunEvals',1200));

%[res,tmp,flag]=fsolve(fitfunc,0.5,optimoptions('fsolve','TolX',1e-8));
% 
% 
% if flag~=1
%     
% 
%     if fitfunc(0.44)<0
%      fitfunc=@(x) -main_func(S,gamma,x(1)*1e-6,B0,TE,theta_B0,rho0,sqrt(dchi_a2/x(1)),R);
%     end
%     
%     res=fminsearch(fitfunc,0.44);
%     
% 
% end


%res=fminsearch(fitfunc,0.8);
% 
% x=linspace(0.01,1);
% for i=1:100
%     y(i)=fitfunc(x(i));
%     
% end
% 
% figure;plot(x,y);
% disp('');
% fitfunc=@(beta,x) main_func(S,gamma,beta(1),B0,TE,theta_B0,rho0,beta(2),R);
% res=nlinfit([],[0,0],fitfunc,[0.50e-6,3]);

function res=main_func(S,gamma,dchi,B0,TE,theta_B0,rho0,a,R)


nTE= length(TE);
res=zeros(1,nTE);
gp=zeros(1,nTE);
p_in=zeros(1,nTE);
    for i=1:length(TE)
        
        p_in(i)=phi_in(gamma,dchi,B0,TE(i),theta_B0);
        
        gp(i)=0.5*gamma*B0*dchi*TE(i)*sin(theta_B0*pi/180)^2/1000;
        
        % only valid for small angle
        % res(i)=real(S(i))*sin(p_in(i))-imag(S(i))*cos(p_in(i))...
        %      -pi*rho0(i)*sin(p_in(i))*(R^2+gp(i)^2*a^4/4/R^2-a^2-gp(i)^2*a^2/4);
        %
        res(i)=real(S(i))*sin(p_in(i))-imag(S(i))*cos(p_in(i))...
            -TotalSignal_Tissue_PerpPlane(a,dchi,B0,TE(i),theta_B0,R)*rho0(i)*sin(p_in(i));
        
    end

return;
res(3)=(rho0(2)*real(S(1))-rho0(1)*real(S(2)))*sin(p_in(1))*sin(p_in(2))...
       -pi*rho0(1)*rho0(2)*sin(p_in(1))*sin(p_in(2))*a^2/4*(1-a^2/R^2)*(gp(2)^2-gp(1)^2)...
       -rho0(2)*imag(S(1))*sin(p_in(2))*cos(p_in(1))...
       +rho0(1)*imag(S(2))*sin(p_in(1))*cos(p_in(2));


function ps=phi_in(gamma,dchi,B0,TE,theta_B0)

g=0.5*gamma*dchi*B0*TE/1000;

theta_B0=theta_B0*pi/180;

    
ps=g/3*(3*cos(theta_B0)^2-1);  %consistent with Siemens; opposite to Hsieh
