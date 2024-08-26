function [pv,v]=PV_V4TOF_PC(TOF,PC, TR,T1f,T1s,thk, FA,VENC)
% [pv,v]=PV_V4TOF_PC(TOF,PC)
% TOF: signal enhancement; ratio compared to the no flow condition
% PC: phase contrast in degrees
% T1d: T1 of the flowing spins in s
% T1s: T1 of the static spins in s
% thk: thickness of the slice in cm
% FA: flip angle in degree
% VENC: venc in cm/s
% pv: partial volume
% v: flow velocity; cm/s
% tof_from_fl_fq: bool. if tof is acquired from the fl_fq scan, then flow
% enhanced signal are acquired once every two TR's.

%% First do some testing
%{
TR = 0.027; % s
T1f = 2.643;
T1s = 2.643;
thk = 0.2; % cm
FA = 13; %
v=0.3; %cm/s
VENC = 1; % cm/s
pv=0.3;

v_meas=V_with_PV(v,thk,TR,FA,T1f,T1s,pv,VENC);
TOF = total_signal(v,thk,TR,FA,T1f,T1s,pv);

%}



TOF=double(TOF);
PC=double(PC);
TR=double(TR);
T1f=double(T1f);
T1s=double(T1s);
thk=double(thk);
FA=double(FA);
VENC=double(VENC);
v_meas=PC*VENC/180;

myfunc=@(x) [V_with_PV(x(2),thk,TR,FA,T1f,T1s,x(1),VENC)-v_meas,total_signal(x(2),thk,TR,FA,T1f,T1s,x(1),1)-TOF]; % x(1) pvf; x(2) v.


res=fsolve(myfunc, [0.42,0.45],optimoptions('fsolve','Display','off','MaxIter',19));           

pv=res(1);
v=res(2);
           
function v_meas=V_with_PV(v,thk,TR,fa,T1art,T1pvs,pvf,venc)

if pvf<0
   % pvf=0;
end
phi=0;
th=v/venc*pi/2;
[ratio,spvs,sart]= total_signal(v,thk,TR,fa,T1art,T1pvs,pvf,1);

atmp=sart*pvf/(spvs*(1-pvf)+sart*pvf);

z1=1-atmp+atmp*exp(1i*(phi+th));
z2=1-atmp+atmp*exp(1i*(phi-th));

dph=phase(z1/z2);
v_meas=dph/pi*venc;






function [ratio,spvs,sart]= total_signal(v,thk,TR,fa,T1art,T1pvs,pvf,sz_init)

if ~exist('sz_init','var')
    sz_init=1;
end

if pvf<0
   % pvf=0;
end

spvs=ssFLASH(fa,TR,T1pvs);  % pvs
sart=0;
sz=sz_init;
irep=1;
while irep<=ceil(thk/v/TR)
    sart_tmp=sz*sin(fa/180*pi);
    sz=sz*cos(fa/180*pi);    
    sz=1+(sz-1)*exp(-TR/T1art);   

    if irep==ceil(thk/v/TR)
      sart=sart+sart_tmp*(thk-v*TR*(irep-1))/thk;
    else
      sart=sart+sart_tmp*TR*v/thk;  
    end    
    
    irep=irep+1;
end
    
ratio=(sart*pvf+spvs*(1-pvf))/spvs;



