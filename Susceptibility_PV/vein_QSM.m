function [res,flag]=vein_QSM(S,rho0,TE,R,B0,theta_B0)
% S: 1*2
% rho0: 1*2
% TE 1*2; in ms
% theta_B0 in degrees
% R: radias in mm
% res: ppm and mm
% the small angle approach in Hsieh_MRI2015b
% due to digitization error of the total signal, this method give poor
% results when theta_B0>
gamma=42.58e6*2*pi; 

fitfunc=@(x) main_func(S,gamma,x(1)*1e-6,B0,TE,theta_B0,rho0,x(2),R);

%[res,tmp,flag]=fsolve(fitfunc,[0.44,0.1],optimoptions('fsolve','TolX',1e-8));
[res,resnorm,tmp2,flag] =  lsqnonlin(fitfunc, [0.44,0.1],[0.2,0],[0.8,1],optimoptions('lsqnonlin','TolFun',1e-12,'MaxFunEvals',1200));


% fitfunc=@(beta,x) main_func(S,gamma,beta(1),B0,TE,theta_B0,rho0,beta(2),R);
% res=nlinfit([],[0,0],fitfunc,[0.50e-6,3]);

function res=main_func(S,gamma,dchi,B0,TE,theta_B0,rho0,a,R)


nTE= length(TE);
res=zeros(1,nTE);
gp=zeros(1,nTE);
p_in=zeros(1,nTE);
for i=1:nTE
    
    p_in(i)=phi_in(gamma,dchi,B0,TE(i),theta_B0);
    
    gp(i)=0.5*gamma*B0*dchi*TE(i)*sin(theta_B0*pi/180)^2/1000;
    
  % only valid for small angle 
  % res(i)=real(S(i))*sin(p_in(i))-imag(S(i))*cos(p_in(i))...
  %      -pi*rho0(i)*sin(p_in(i))*(R^2+gp(i)^2*a^4/4/R^2-a^2-gp(i)^2*a^2/4);
  %      
  %  s_tissue=TotalSignal_Tissue_PerpPlane(a,dchi,B0,TE,theta,R)
  res(i)=real(S(i))*sin(p_in(i))-imag(S(i))*cos(p_in(i))...
        -TotalSignal_Tissue_PerpPlane(a,dchi,B0,TE(i),theta_B0,R)*sin(p_in(i))*rho0(i);
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
