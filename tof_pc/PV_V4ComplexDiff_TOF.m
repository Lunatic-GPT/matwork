function [pvf,v,ok,resnorm]=PV_V4ComplexDiff_TOF(z, v_calc,sart_calc,v_calc_pc,sart_calc_pc,s_static,VENC,lambda,flow_pattern,twoPool)

% z: [deltaZ,z1], where
% deltaZ: complex difference signal normalized by the signal that would be measured if fully occupied by the flowing spin but with no flow enhancement.  
% z1: signal intensity without flow encoding and normalized by the same way

% signals in sart_calc, sart_calc_pc, and s_static are normalized in the same way

% x(1); partial volume fraction
% x(2): mean velocity

% If not two pool (CD alone), then the solution can be obtained when |z(1)|<=|z(2)| and real(z(1))<0;
% 
% 8/15/2017: 
% changed from 2 unknown fitting (v,pv) to 3 unknown fitting (v,pv,center of circle on real axis) 
% for 2 unknown fitting, center of circle was calculated in advance.
% However, no solution for center exists when (real(z(1))>=0 || abs(z(1))>abs(z(2)))
% for now only consider ~twoPool model
if twoPool
    error('Not implemented yet');
end

z1=z(2);  % signal without flow encoding gradient

z2=z(1)+z1; % signal with flow encoding gradient

myfunc=@(x) FitFunc(z1,z2,x(3),x(2),x(1),VENC,v_calc,sart_calc,v_calc_pc,sart_calc_pc,flow_pattern);

%[res,resnorm,tmp2,ok]=lsqnonlin(myfunc, [0.1,VENC/2,0.5*z1],[0,0,0],[2,VENC*2,z1],optimoptions('lsqnonlin','TolFun',1e-12,'MaxFunctionEvaluations',1200));
[res,resnorm,tmp2,ok]=lsqnonlin(myfunc, [0.1,VENC/2,0.5*z1],[0,0,0],[2,VENC*2,z1],optimoptions('lsqnonlin','TolFun',1e-12,'MaxFunEvals',1200));
%MaxFunEvals
 
pvf=res(1);
v=res(2);
%disp(myfunc(res));

    
% code below for old implementation
% 
% if ~twoPool && (real(z(1))>=0 || abs(z(1))>abs(z(2)))
%  
%     
%     disp('PV_V4ComplexDiff_TOF: Can not solve');
%     pvf=0;
%     v=0;
%     ok=false;
%     return;
% end
% 
% 
% if ~twoPool
% myfunc=@(x) [CmplxDiff(x(2),x(1),VENC,v_calc_pc,sart_calc_pc,flow_pattern,1)-real(z(1)),...
%              CmplxDiff(x(2),x(1),VENC,v_calc_pc,sart_calc_pc,flow_pattern,2)-imag(z(1))]; % x(1) pvf; x(2) v.
%          
% else
%     
% myfunc=@(x) [CmplxDiff(x(2),x(1),VENC,v_calc_pc,sart_calc_pc,flow_pattern,1)-real(z(1)),...
%              CmplxDiff(x(2),x(1),VENC,v_calc_pc,sart_calc_pc,flow_pattern,2)-imag(z(1)),...
%              (M_noFlowEncoding(x(2),x(1),lambda,v_calc,sart_calc,s_static)-z(2))*2]; % x(1) pvf; x(2) v.
%          
% end
%          
% 
% %%
% ang1=pi-angle0_2pi(z(1));
% 
% func_tmp=@(x) (angle0_2pi(z(1)+z(2)-x)-pi+2*ang1);
% 
% xint=fsolve(func_tmp,0.5*z(2));   % the center of circle on x-axis;
% 
% v0=angle0_2pi(z(1)+z(2)-xint)/pi*VENC;
% 
% 
% s1=M_noFlowEncoding(v0,1,lambda,v_calc,sart_calc,s_static(1));
% s2=M_noFlowEncoding(v0,0,lambda,v_calc,sart_calc,s_static(1));
% 
% r=xint/z(2);
% 
% pvf0=s2*(1-r)/(r*s1-r*s2+s2);
% 
% if pvf0>1  || pvf0<0 || v0<0
%     res=[0,0];   
%     resid=myfunc(res)/z(2);
%     ok=false;
%     
% else
%     
%     %%
%     
%     %disp(myfunc([pvf0,v0,scale0]));
%     %try
%     %res=lsqnonlin(myfunc, [pvf0,v0],[0,0],[1,5],optimoptions('lsqnonlin','FunctionTolerance',1e-12));
%     res=lsqnonlin(myfunc, [pvf0,v0],[0,0],[1,10],optimoptions('lsqnonlin','TolFun',1e-12));
%     %
%     % catch
%     %     res=[0,0,0];
%     % end
%     
%     resid=myfunc(res)/z(2);
%     if any(abs(resid)>0.03)
%         ok=false;
%     else
%         ok=true;
%     end
%     
% end
% %%
% 
% if exist('verbose','var') && verbose
%     fprintf('Residual =');
%     fprintf(' %f; ',resid);
%     fprintf('\n');
% end
% 
% pvf=res(1);
% v=res(2);

%%
 
function res=FitFunc(z1,z2,c,v,pv,VENC,v_calc,sart_calc,v_calc_pc,sart_calc_pc,flow_pattern)

sflow = FlowSignal(v,pv,VENC,v_calc_pc,sart_calc_pc,flow_pattern);

s=interp1(v_calc,sart_calc,v)*pv;

res(1)=real(z2-c-sflow);
res(2)=imag(z2-c-sflow);
res(3)=(c+s-z1);

%res = abs(z2-c+sflow).^2+(c+s-z1)^2;



function res = M_noFlowEncoding(v,pvf,lambda,v_calc,sart_calc,s_static)

sart= interp1(v_calc,sart_calc,v);

res=(1-pvf)*s_static/lambda+pvf*sart;


function res=CmplxDiff(v,pvf,venc,v_calc,sart_calc,flow_pattern,re_img)

% lambda: the spin density ratio between flowing and static spins 
% maximum velocity is twice the mean velocity for laminar flow

if strcmp(flow_pattern,'laminar')
r=linspace(0,1,1000);
vloc=v*2*(1-r.^2);  %maximum velocity is twice the mean velocity
elseif strcmp(flow_pattern,'plug')
    r=1;
   vloc=v;
else
    error('Unknown pattern');
end

th=vloc/venc*pi;

sart= interp1(v_calc,sart_calc,vloc);

%sart= TOF_signalIntensity(v,thk,TR,fa,T1(1),flow_pattern,pulse_profile);
%s_static=TOF_signalIntensity(0,thk,TR,fa,T1(2),flow_pattern,pulse_profile);

cd=pvf*sum(sart.*(exp(1i*th)-1).*r)/sum(r);


if re_img==1
    res=real(cd);
elseif re_img==2
    res=imag(cd);
else 
    res=cd;
end


function res=FlowSignal(v,pvf,venc,v_calc,sart_calc,flow_pattern)

% lambda: the spin density ratio between flowing and static spins 
% maximum velocity is twice the mean velocity for laminar flow

if strcmp(flow_pattern,'laminar')
r=linspace(0,1,1000);
vloc=v*2*(1-r.^2);  %maximum velocity is twice the mean velocity
elseif strcmp(flow_pattern,'plug')
    r=1;
   vloc=v;
else
    error('Unknown pattern');
end

th=vloc/venc*pi;

sart= interp1(v_calc,sart_calc,vloc);

%sart= TOF_signalIntensity(v,thk,TR,fa,T1(1),flow_pattern,pulse_profile);
%s_static=TOF_signalIntensity(0,thk,TR,fa,T1(2),flow_pattern,pulse_profile);

res=pvf*sum(sart.*exp(1i*th).*r)/sum(r);




