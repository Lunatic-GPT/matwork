function [TOF,PC]=TOF_PC4PV_V(pv,v, TR,T1f,T1s,thk, FA,VENC)
% [TOF,PC]=TOF_PC4PV_V(pv,v, TR,T1f,T1s,thk, FA,VENC)
% TOF: signal enhancement; ratio compared to the no flow condition
% PC: phase contrast in degrees
% T1d: T1 of the flowing spins in s
% T1s: T1 of the static spins in s
% thk: thickness of the slice in cm
% FA: flip angle in degree
% VENC: venc in cm/s
% pv: partial volume
% v: flow velocity; cm/s
% 

%% First do some testing

v_meas=V_with_PV(v,thk,TR,FA,T1f,T1s,pv,VENC);
TOF = total_signal(v,thk,TR,FA,T1f,T1s,pv);

PC=v_meas/VENC*180;

           
function v_meas=V_with_PV(v,thk,TR,fa,T1art,T1pvs,pvf,venc)

phi=0;
th=v/venc*pi/2;
[ratio,spvs,sart]= total_signal(v,thk,TR,fa,T1art,T1pvs,pvf);

atmp=sart*pvf/(spvs*(1-pvf)+sart*pvf);

z1=1-atmp+atmp*exp(1i*(phi+th));
z2=1-atmp+atmp*exp(1i*(phi-th));

dph=phase(z1/z2);
v_meas=dph/pi*venc;






function [ratio,spvs,sart]= total_signal(v,thk,TR,fa,T1art,T1pvs,pvf,sz_init)

if ~exist('sz_init','var')
    sz_init=1;
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



