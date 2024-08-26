function [pvf,v,ok]=PV_V4PC_Hamilton(z, v_calc,sart_calc,v_calc_pc,sart_calc_pc,s_static,VENC,lambda,flow_pattern,verbose)

% z: [z1,z2] or [z1,z2,zref], where
% z1: complex signal without phase encoding gradient
% z2: complex signal with phase encoding gradient
% zref: complex signal without phase encoding gradient in voxels containing
% only flowing spin with v=0;
% s_static: signal from static tissue without partial volume effect

% x(1); partial volume fraction
% x(2): mean velocity
% x(3): scaling to match signal


ph=angle(z(1));
m1=abs(z(1));
m2=z(2)*exp(-1i*ph);


if length(z)==2
myfunc=@(x) [M_noFlowEncoding(x(2),x(1),x(3),lambda,v_calc,sart_calc,s_static)-m1, ...
             M_FlowEncoding(x(2),x(1),x(3),VENC,lambda,v_calc_pc,sart_calc_pc,s_static,flow_pattern,1)-real(m2),...
             M_FlowEncoding(x(2),x(1),x(3),VENC,lambda,v_calc_pc,sart_calc_pc,s_static,flow_pattern,2)-imag(m2),...
             ]; % x(1) pvf; x(2) v.
         
else
    
    mref=abs(z(3));
    
    myfunc=@(x) [M_noFlowEncoding(x(2),x(1),x(3),lambda,v_calc,sart_calc,s_static)-m1, ...
             M_FlowEncoding(x(2),x(1),x(3),VENC,lambda,v_calc_pc,sart_calc_pc,s_static,flow_pattern,1)-real(m2),...
             M_FlowEncoding(x(2),x(1),x(3),VENC,lambda,v_calc_pc,sart_calc_pc,s_static,flow_pattern,2)-imag(m2),...
             M_noFlowEncoding(0,1,x(3),lambda,v_calc,sart_calc,s_static)-mref
             ]; % x(1) pvf; x(2) v.
end
         

%%
ang1=pi-angle((z(2)-z(1)));

func_tmp=@(x) (angle(z(2)-x)-pi+2*ang1);

xint=fsolve(func_tmp,0.5);

v0=angle(z(2)-xint)/pi*VENC;

pvf0=1-xint/m1;

s1=M_noFlowEncoding(v0,1,1,lambda,v_calc,sart_calc,s_static);
s2=M_noFlowEncoding(v0,0,1,lambda,v_calc,sart_calc,s_static);

r=xint/m1;

pvf0=s2*(1-r)/(r*s1-r*s2+s2);



if pvf0>1  || pvf0<0 || v0<0
    res=[0,0,0];
else

%%d
scale0=m1/M_noFlowEncoding(v0,pvf0,1,lambda,v_calc,sart_calc,s_static);

%disp(myfunc([pvf0,v0,scale0]));
%try
%res=lsqnonlin(myfunc, [pvf0,v0,scale0],[0,0,0],[1,10,scale0*10],optimoptions('lsqnonlin','StepTolerance',1e-6,'FunctionTolerance',1e-12)); 
res=lsqnonlin(myfunc, [pvf0,v0,scale0],[0,0,0],[1,10,scale0*10],optimoptions('lsqnonlin','TolFun',1e-12)); 
% 
% catch
%     res=[0,0,0];
% end

end
%%

resid=myfunc(res)/m1;
if any(abs(resid)>0.03) 
    ok=false;
else
    ok=true;
end

if exist('verbose','var') && verbose
    fprintf('Residual =');
    fprintf(' %f; ',resid);
    fprintf('\n');
end

pvf=res(1);
v=res(2);

%%
 

function res = M_noFlowEncoding(v,pvf,scale,lambda,v_calc,sart_calc,s_static)

sart= interp1(v_calc,sart_calc,v);

res=(1-pvf)*s_static+pvf*lambda*sart;

res=res*scale;

function res=M_FlowEncoding(v,pvf,scale,venc,lambda,v_calc,sart_calc,s_static,flow_pattern,re_img)

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

z1=(1-pvf)*s_static+pvf*lambda*sum(sart.*exp(1i*th).*r)/sum(r);

z1=scale*z1;

if re_img==1
    res=real(z1);
elseif re_img==2
    res=imag(z1);
end




function sart_total= TOF_signalIntensityFast(v,v_calc,sart_calc,flow_pattern)
%TOF_signalIntensity(v,thk,TR,fa,T1,flow_pattern,pulse_profile,sz_init)
% flow_pattern: laminar or plug
% pulse_profile: square or sinc
% thk in cm
% v: mean velocity in cm/s
% TR: in s
% T1: in s
% sz_init: default 1

if ~exist('sz_init','var')
    sz_init=1;
end

if strcmp(flow_pattern,'laminar')
r=linspace(0,1,100);
vloc=v*2*(1-r.^2);  %maximum velocity is twice the mean velocity
elseif strcmp(flow_pattern,'plug')
    r=1;
   vloc=v;
else
    error('Unknown pattern');
end


tic;
sart=zeros(1,length(r));
for i=1:length(r)
  
 %  sart(i)=total_signal(vloc(i), TR,pos,fa,T1,sz_init);
  sart(i)=interp1(v_calc,sart_calc,vloc(i));

end

sart_total=sum(sart.*r)/sum(r);  %weighted average




function sart=total_signal(v,TR,thk,fa,T1,pulse_profile,sz_init)
% pos and fa gives the RF slice selection profile;
% for an ideal boxcar profile use pos=[-thk/2 thk/2], fa=[fa,fa];


if ~exist('sz_init','var')
    sz_init=1;
end

if strcmp(pulse_profile,'square')
    pos=[-thk/2 thk/2];
    fa=[fa,fa];
elseif strcmp(pulse_profile,'sinc')
   
 %   pos=ri('SINC_Profile_Normalized4FA25.mat','','','d');
 pos= linspace(-6,6,400);
 faNorm=[3.0e-05 3.0e-05 2.9e-05 2.5e-05 2.0e-05 1.4e-05 5.9e-06 2.3e-06 1.1e-05 1.9e-05 2.6e-05 3.2e-05 3.6e-05 3.8e-05 3.8e-05 3.5e-05 3.0e-05 2.3e-05 1.5e-05 4.8e-06 5.7e-06 1.6e-05 2.6e-05 3.5e-05 4.2e-05 4.7e-05 4.9e-05 4.8e-05 4.4e-05 3.7e-05 2.7e-05 1.6e-05 2.6e-06 1.1e-05 2.5e-05 3.7e-05 4.8e-05 5.7e-05 6.2e-05 6.3e-05 6.1e-05 5.5e-05 4.5e-05 3.2e-05 1.6e-05 1.7e-06 2.0e-05 3.8e-05 5.4e-05 6.8e-05 7.8e-05 8.4e-05 8.5e-05 8.0e-05 7.1e-05 5.6e-05 3.7e-05 1.5e-05 9.8e-06 3.5e-05 5.9e-05 8.1e-05 9.9e-05 1.1e-04 1.2e-04 1.2e-04 1.1e-04 9.3e-05 7.1e-05 4.3e-05 1.0e-05 2.5e-05 6.1e-05 9.5e-05 1.2e-04 1.5e-04 1.6e-04 1.7e-04 1.7e-04 1.5e-04 1.3e-04 9.0e-05 4.6e-05 4.3e-06 5.8e-05 1.1e-04 1.6e-04 2.1e-04 2.4e-04 2.6e-04 2.7e-04 2.6e-04 2.3e-04 1.9e-04 1.3e-04 5.5e-05 2.6e-05 1.1e-04 2.0e-04 2.8e-04 3.4e-04 3.9e-04 4.2e-04 4.2e-04 3.9e-04 3.3e-04 2.4e-04 1.2e-04 2.3e-05 1.8e-04 3.5e-04 5.1e-04 6.6e-04 7.9e-04 8.8e-04 9.3e-04 9.3e-04 8.7e-04 7.5e-04 5.7e-04 3.3e-04 5.0e-05 2.7e-04 6.0e-04 9.2e-04 1.2e-03 1.4e-03 1.6e-03 1.6e-03 1.5e-03 1.3e-03 8.7e-04 2.7e-04 4.9e-04 1.4e-03 2.4e-03 3.6e-03 4.7e-03 5.7e-03 6.6e-03 7.2e-03 7.5e-03 7.1e-03 6.1e-03 4.1e-03 1.2e-03 3.0e-03 8.6e-03 1.6e-02 2.5e-02 3.5e-02 4.8e-02 6.3e-02 8.0e-02 1.0e-01 1.2e-01 1.5e-01 1.7e-01 2.0e-01 2.3e-01 2.6e-01 3.0e-01 3.3e-01 3.7e-01 4.1e-01 4.4e-01 4.8e-01 5.2e-01 5.6e-01 6.0e-01 6.4e-01 6.7e-01 7.1e-01 7.4e-01 7.7e-01 8.0e-01 8.3e-01 8.6e-01 8.8e-01 9.0e-01 9.2e-01 9.4e-01 9.5e-01 9.7e-01 9.8e-01 9.8e-01 9.9e-01 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 9.9e-01 9.8e-01 9.8e-01 9.7e-01 9.5e-01 9.4e-01 9.2e-01 9.0e-01 8.8e-01 8.6e-01 8.3e-01 8.0e-01 7.7e-01 7.4e-01 7.1e-01 6.7e-01 6.4e-01 6.0e-01 5.6e-01 5.2e-01 4.8e-01 4.4e-01 4.1e-01 3.7e-01 3.3e-01 3.0e-01 2.6e-01 2.3e-01 2.0e-01 1.7e-01 1.5e-01 1.2e-01 1.0e-01 8.0e-02 6.3e-02 4.8e-02 3.5e-02 2.5e-02 1.6e-02 8.6e-03 3.0e-03 1.2e-03 4.1e-03 6.1e-03 7.1e-03 7.5e-03 7.2e-03 6.6e-03 5.7e-03 4.7e-03 3.6e-03 2.4e-03 1.4e-03 4.9e-04 2.7e-04 8.7e-04 1.3e-03 1.5e-03 1.6e-03 1.6e-03 1.4e-03 1.2e-03 9.2e-04 6.0e-04 2.7e-04 5.0e-05 3.3e-04 5.7e-04 7.5e-04 8.7e-04 9.3e-04 9.3e-04 8.8e-04 7.9e-04 6.6e-04 5.1e-04 3.5e-04 1.8e-04 2.3e-05 1.2e-04 2.4e-04 3.3e-04 3.9e-04 4.2e-04 4.2e-04 3.9e-04 3.4e-04 2.8e-04 2.0e-04 1.1e-04 2.6e-05 5.5e-05 1.3e-04 1.9e-04 2.3e-04 2.6e-04 2.7e-04 2.6e-04 2.4e-04 2.1e-04 1.6e-04 1.1e-04 5.8e-05 4.3e-06 4.6e-05 9.0e-05 1.3e-04 1.5e-04 1.7e-04 1.7e-04 1.6e-04 1.5e-04 1.2e-04 9.5e-05 6.1e-05 2.5e-05 1.0e-05 4.3e-05 7.1e-05 9.3e-05 1.1e-04 1.2e-04 1.2e-04 1.1e-04 9.9e-05 8.1e-05 5.9e-05 3.5e-05 9.8e-06 1.5e-05 3.7e-05 5.6e-05 7.1e-05 8.0e-05 8.5e-05 8.4e-05 7.8e-05 6.8e-05 5.4e-05 3.8e-05 2.0e-05 1.7e-06 1.6e-05 3.2e-05 4.5e-05 5.5e-05 6.1e-05 6.3e-05 6.2e-05 5.7e-05 4.8e-05 3.7e-05 2.5e-05 1.1e-05 2.6e-06 1.6e-05 2.7e-05 3.7e-05 4.4e-05 4.8e-05 4.9e-05 4.7e-05 4.2e-05 3.5e-05 2.6e-05 1.6e-05 5.7e-06 4.8e-06 1.5e-05 2.3e-05 3.0e-05 3.5e-05 3.8e-05 3.8e-05 3.6e-05 3.2e-05 2.6e-05 1.9e-05 1.1e-05 2.3e-06 5.9e-06 1.4e-05 2.0e-05 2.5e-05 2.9e-05 3.0e-05 3.0e-05];
 
 pos=pos(101:300);
 faNorm=faNorm(101:300);
 %faNorm=ri('SINC_Profile_Normalized4FA25.mat','','','angle');
    
    pos=pos*thk/2;
    fa=fa*faNorm;
    
else
    error('Unknown pulse profile');
end
%T1art=0.2;
%T1tissue=1.4;
pos=pos-min(pos);  % start from 0
thk2=max(pos);

sz=sz_init;

if v==0
    
    fa=fa(1:4:end);
    
    for i=1:length(fa)
     sart_tmp(i)=ssFLASH(fa(i),TR,T1,0,1);
    end
    sart=mean(sart_tmp);
    
sart=sart*thk2/thk;
    return;
end

TotalRF =ceil(thk2/v/TR);
sart_tmp=zeros(1,TotalRF);

faa=zeros(1,TotalRF);

for irep=1:TotalRF
   
    if irep<TotalRF
        z=(irep-0.5)*v*TR;
    else
        z=(thk2+(irep-1)*v*TR)/2;
    end
    faa(irep)=interp1(pos,fa,z);
    
  % faa(irep)=fa(1)+(fa(2)-fa(1))/thk2*z;
   
    
end

for irep=1:TotalRF
    sart_tmp(irep)=sz*sin(faa(irep)/180*pi);
    sz=sz*cos(faa(irep)/180*pi);    
    sz=1+(sz-1)*exp(-TR/T1);   
end


sart=0;
for irep=1:TotalRF
    
    if irep==ceil(thk2/v/TR)
      sart=sart+sart_tmp(irep)*(thk2-v*TR*(irep-1))/thk2;
    else
      sart=sart+sart_tmp(irep)*TR*v/thk2;  
    end
    
end
sart=sart*thk2/thk;

