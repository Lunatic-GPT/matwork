function res=vein_QSM_largeAngle(S,rho0,TE,R,B0,theta_B0)
% S: 1*2
% rho0: 1*2
% TE 1*2;

gamma=42.58e6*2*pi; 



%fitfunc=@(x) main_func(S,gamma,x(1),B0,TE,theta_B0,rho0,1.83,R);

%res=fsolve(fitfunc,[0.50e-6]);

fitfunc=@(beta,x) main_func(S,gamma,beta(1),B0,TE,theta_B0,rho0,1.83,R);
res=nlinfit([],[0,0],fitfunc,[0.50e-6]);

function res=main_func(S,gamma,dchi,B0,TE,theta_B0,rho0,a,R)


res=zeros(1,3);
gp=zeros(1,2);
p_in=zeros(1,2);
for i=1:2
    
    p_in(i)=phi_in(gamma,dchi,B0,TE(i),theta_B0);
    
    gp(i)=0.5*gamma*B0*dchi*TE(i)*sin(theta_B0)^2;
    
    res(i)=real(S(i))*sin(p_in(i))-imag(S(1))*cos(p_in(i))...
        -pi*rho0(i)*sin(p_in(i))*(R^2-gp(i)^2*a^4/4/R^2-a^2-gp(i)^2*a^2/4);
    
end

res=res(1:2);
return;
res(3)=(rho0(2)*real(S(1))-rho0(1)*real(S(2)))*sin(p_in(1))*sin(p_in(2))...
       -pi*rho0(1)*rho0(2)*sin(p_in(1))*sin(p_in(2))*a^2/4*(1-a^2/R^2)*(gp(2)^2-gp(1)^2)...
       -rho0(2)*imag(S(1))*sin(p_in(2))*cos(p_in(1))...
       +rho0(1)*imag(S(2))*sin(p_in(1))*cos(p_in(2));


function ps=phi_in(gamma,dchi,B0,TE,theta_B0)

g=0.5*gamma*dchi*B0*TE/1000;

theta_B0=theta_B0*pi/180;

    
ps=-g/3*(3*cos(theta_B0)^2-1);
