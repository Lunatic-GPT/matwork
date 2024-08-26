function v_meas=V_with_PV(v,thk,TR,fa,T1,TE,T2,pvf,venc,lambda,flow_pattern,pulse_profile)

% lambda: the spin density ratio between flowing and static spins 
if length(T1)==1
    T1 = T1*[1,1];
end

if length(T2)==1
    T2 = T2*[1,1];
end
 %maximum velocity is twice the mean velocity

 if strcmp(flow_pattern,'laminar')
r=linspace(0,1,100);
vloc=v*2*(1-r.^2);  %maximum velocity is twice the mean velocity
elseif strcmp(flow_pattern,'plug')
    r=1;
   vloc=v;
else
    error('Unknown pattern');
end

th=vloc/venc*pi/2;


sart= TOF_signalIntensity(v,thk,TR,fa,T1(1),flow_pattern,pulse_profile);
spvs=TOF_signalIntensity(0,thk,TR,fa,T1(2),flow_pattern,pulse_profile);



z1=(1-pvf)*spvs*exp(-TE/T2(2))+pvf*lambda*sum(sart.*exp(1i*th).*r)/sum(r)*exp(-TE/T2(1));
z2=(1-pvf)*spvs*exp(-TE/T2(2))+pvf*lambda*sum(sart.*exp(-1i*th).*r)/sum(r)*exp(-TE/T2(1));


dph=phase(z1/z2);
v_meas=dph/pi*venc;





