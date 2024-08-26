function [v_calc,sart_calc,v_calc_pc,sart_calc_pc,s_static]=Prep_PV_V4TOF2Sur_PC(TR,T1,TE,T2,thk,FA,flow_pattern,pulse_profile,exp_sel)
%  [v_calc,sart_calc,v_calc_pc,sart_calc_pc,s_static]=Prep_PV_V4TOF2Sur_PC(TR,T1,TE,T2,thk,FA,flow_pattern,pulse_profile,exp_sel)
% calculate the interpolation curves for PV_V4TOF2Sur_PC calculatioin.
% TR: 1*2 vector; TR of the TOF and PC scans
% T1: 1*2 vector; T1 of the flowing and static spins  in s
% TE: 1*2 vector; TE of the TOF and PC scans
% T2: 1*2 vector; T2 of the flowing and static spins in s
% FA: 1*2 vector; flip angle of the TOF and PC scans in degrees
% thk: thickness of the slice in cm
% flow_pattern: laminar or plut
% pulse_profile: sinc or square
% exp_sel: fl_fq_retroz_mb or fl_fq_retroz
% s_static : a 1*2 array containing signals for the static spins without
% partial volume effects in the TOF and PC scans; The results assumes the
% same water density for static and flowing spins.
% tilt of the vessel relative to the normal direction (in degree); default 0
if ~exist('exp_sel','var')
    exp_sel='fl_fq_retroZ_mb';
end

%%

fname=param2name(TR,T1,TE,T2,thk,FA,flow_pattern,pulse_profile,exp_sel);

if exist(fullfile(toml,'TOF_PC','Prep_PV_V4TOF2Sur_PC_Save',fname),'file')
    load(fullfile(toml,'TOF_PC','Prep_PV_V4TOF2Sur_PC_Save',fname));
    
%     v_calc=tmp.v_calc;
%     sart_calc=tmp.sart_calc;
%     v_calc_pc=tmp.v_calc_pc;
%     sart_calc_pc=tmp.sart_calc_pc;
%     s_static=tmp.s_static;
%     
    return;
end

same_pc_tof=0;
if length(TR)==1
    same_pc_tof=1;
    TR = TR*[1,1];
end

if length(TE)==1
    TE = TE*[1,1];
end

if length(T1)==1
    T1 = T1*[1,1];
end

if length(T2)==1
    T2 = T2*[1,1];
end

if length(FA)==1
    FA = FA*[1,1];
end

TR=double(TR);
T1=double(T1);
T2=double(T2);
thk=double(thk);
FA=double(FA);


%%

v_calc_pre=[0:0.01:0.05,0.06:0.02:10,11:90];
%v_calc_pre=0.7692;
%v_calc_pre=60;
disp('precalculate the v vs signal curves');
sart_calc_pre=zeros(length(v_calc_pre),2);

for i=1:length(v_calc_pre)
    if same_pc_tof
        sart_calc_pre(i,:)=total_signal_new(v_calc_pre(i),TR(1),thk,FA(1),T1(1),pulse_profile,1,exp_sel);
    else
        for scan=1:2        
            sart_calc_pre(i,scan)=total_signal_new(v_calc_pre(i),TR(scan),thk,FA(scan),T1(1),pulse_profile,1,exp_sel);
        end
    end
end
% TOF_signalIntensity(v,v_calc,sart_calc,flow_pattern)
v_calc=[0:0.1:15,16:90];
%v_calc=v_calc_pre;
for i=1:length(v_calc)
    sart_calc(i)= TOF_signalIntensityFast(v_calc(i),v_calc_pre,sart_calc_pre(:,1),flow_pattern)*exp(-TE(1)/T2(1));
end

for scan=1:2
    s_static(scan)=total_signal_new(0,TR(scan),thk,FA(scan),T1(2),pulse_profile,1,exp_sel)*exp(-TE(scan)/T2(2));  % was using total_signal; should be the same as total_signal
end


sart_calc_pc=sart_calc_pre(:,2)*exp(-TE(2)/T2(1));
v_calc_pc=v_calc_pre;


s0=sart_calc_pc(1);
sart_calc_pc=sart_calc_pc/s0;
s_static(2)=s_static(2)/s0;

s0=sart_calc(1);
sart_calc=sart_calc/s0;
s_static(1)=s_static(1)/s0;

if same_pc_tof
    s_static = s_static(1);
end


save(fullfile(toml,'TOF_PC','Prep_PV_V4TOF2Sur_PC_Save',fname),'v_calc','sart_calc','v_calc_pc','sart_calc_pc','s_static');


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

function sart=total_signal_new(v,TR,thk,fa,T1,pulse_profile,sz_init,exp_sel)
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
    if strcmp(exp_sel,'fl_fq_retroZ_mb')
        pos= linspace(-6,6,400);
        faNorm=[3.0e-05 3.0e-05 2.9e-05 2.5e-05 2.0e-05 1.4e-05 5.9e-06 2.3e-06 1.1e-05 1.9e-05 2.6e-05 3.2e-05 3.6e-05 3.8e-05 3.8e-05 3.5e-05 3.0e-05 2.3e-05 1.5e-05 4.8e-06 5.7e-06 1.6e-05 2.6e-05 3.5e-05 4.2e-05 4.7e-05 4.9e-05 4.8e-05 4.4e-05 3.7e-05 2.7e-05 1.6e-05 2.6e-06 1.1e-05 2.5e-05 3.7e-05 4.8e-05 5.7e-05 6.2e-05 6.3e-05 6.1e-05 5.5e-05 4.5e-05 3.2e-05 1.6e-05 1.7e-06 2.0e-05 3.8e-05 5.4e-05 6.8e-05 7.8e-05 8.4e-05 8.5e-05 8.0e-05 7.1e-05 5.6e-05 3.7e-05 1.5e-05 9.8e-06 3.5e-05 5.9e-05 8.1e-05 9.9e-05 1.1e-04 1.2e-04 1.2e-04 1.1e-04 9.3e-05 7.1e-05 4.3e-05 1.0e-05 2.5e-05 6.1e-05 9.5e-05 1.2e-04 1.5e-04 1.6e-04 1.7e-04 1.7e-04 1.5e-04 1.3e-04 9.0e-05 4.6e-05 4.3e-06 5.8e-05 1.1e-04 1.6e-04 2.1e-04 2.4e-04 2.6e-04 2.7e-04 2.6e-04 2.3e-04 1.9e-04 1.3e-04 5.5e-05 2.6e-05 1.1e-04 2.0e-04 2.8e-04 3.4e-04 3.9e-04 4.2e-04 4.2e-04 3.9e-04 3.3e-04 2.4e-04 1.2e-04 2.3e-05 1.8e-04 3.5e-04 5.1e-04 6.6e-04 7.9e-04 8.8e-04 9.3e-04 9.3e-04 8.7e-04 7.5e-04 5.7e-04 3.3e-04 5.0e-05 2.7e-04 6.0e-04 9.2e-04 1.2e-03 1.4e-03 1.6e-03 1.6e-03 1.5e-03 1.3e-03 8.7e-04 2.7e-04 4.9e-04 1.4e-03 2.4e-03 3.6e-03 4.7e-03 5.7e-03 6.6e-03 7.2e-03 7.5e-03 7.1e-03 6.1e-03 4.1e-03 1.2e-03 3.0e-03 8.6e-03 1.6e-02 2.5e-02 3.5e-02 4.8e-02 6.3e-02 8.0e-02 1.0e-01 1.2e-01 1.5e-01 1.7e-01 2.0e-01 2.3e-01 2.6e-01 3.0e-01 3.3e-01 3.7e-01 4.1e-01 4.4e-01 4.8e-01 5.2e-01 5.6e-01 6.0e-01 6.4e-01 6.7e-01 7.1e-01 7.4e-01 7.7e-01 8.0e-01 8.3e-01 8.6e-01 8.8e-01 9.0e-01 9.2e-01 9.4e-01 9.5e-01 9.7e-01 9.8e-01 9.8e-01 9.9e-01 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 9.9e-01 9.8e-01 9.8e-01 9.7e-01 9.5e-01 9.4e-01 9.2e-01 9.0e-01 8.8e-01 8.6e-01 8.3e-01 8.0e-01 7.7e-01 7.4e-01 7.1e-01 6.7e-01 6.4e-01 6.0e-01 5.6e-01 5.2e-01 4.8e-01 4.4e-01 4.1e-01 3.7e-01 3.3e-01 3.0e-01 2.6e-01 2.3e-01 2.0e-01 1.7e-01 1.5e-01 1.2e-01 1.0e-01 8.0e-02 6.3e-02 4.8e-02 3.5e-02 2.5e-02 1.6e-02 8.6e-03 3.0e-03 1.2e-03 4.1e-03 6.1e-03 7.1e-03 7.5e-03 7.2e-03 6.6e-03 5.7e-03 4.7e-03 3.6e-03 2.4e-03 1.4e-03 4.9e-04 2.7e-04 8.7e-04 1.3e-03 1.5e-03 1.6e-03 1.6e-03 1.4e-03 1.2e-03 9.2e-04 6.0e-04 2.7e-04 5.0e-05 3.3e-04 5.7e-04 7.5e-04 8.7e-04 9.3e-04 9.3e-04 8.8e-04 7.9e-04 6.6e-04 5.1e-04 3.5e-04 1.8e-04 2.3e-05 1.2e-04 2.4e-04 3.3e-04 3.9e-04 4.2e-04 4.2e-04 3.9e-04 3.4e-04 2.8e-04 2.0e-04 1.1e-04 2.6e-05 5.5e-05 1.3e-04 1.9e-04 2.3e-04 2.6e-04 2.7e-04 2.6e-04 2.4e-04 2.1e-04 1.6e-04 1.1e-04 5.8e-05 4.3e-06 4.6e-05 9.0e-05 1.3e-04 1.5e-04 1.7e-04 1.7e-04 1.6e-04 1.5e-04 1.2e-04 9.5e-05 6.1e-05 2.5e-05 1.0e-05 4.3e-05 7.1e-05 9.3e-05 1.1e-04 1.2e-04 1.2e-04 1.1e-04 9.9e-05 8.1e-05 5.9e-05 3.5e-05 9.8e-06 1.5e-05 3.7e-05 5.6e-05 7.1e-05 8.0e-05 8.5e-05 8.4e-05 7.8e-05 6.8e-05 5.4e-05 3.8e-05 2.0e-05 1.7e-06 1.6e-05 3.2e-05 4.5e-05 5.5e-05 6.1e-05 6.3e-05 6.2e-05 5.7e-05 4.8e-05 3.7e-05 2.5e-05 1.1e-05 2.6e-06 1.6e-05 2.7e-05 3.7e-05 4.4e-05 4.8e-05 4.9e-05 4.7e-05 4.2e-05 3.5e-05 2.6e-05 1.6e-05 5.7e-06 4.8e-06 1.5e-05 2.3e-05 3.0e-05 3.5e-05 3.8e-05 3.8e-05 3.6e-05 3.2e-05 2.6e-05 1.9e-05 1.1e-05 2.3e-06 5.9e-06 1.4e-05 2.0e-05 2.5e-05 2.9e-05 3.0e-05 3.0e-05];
        pos=pos(101:300);
        faNorm=faNorm(101:300);
        
    elseif strcmp(exp_sel,'fl_fq_retroZ')
        pos=linspace(-4,4,400);
        
        faNorm=[1.725e-03,1.719e-03,1.705e-03,1.681e-03,1.648e-03,1.605e-03,1.552e-03,1.488e-03,1.414e-03,1.329e-03,1.233e-03,1.126e-03,1.007e-03,8.766e-04,7.347e-04,5.812e-04,4.162e-04,2.398e-04,5.204e-05,1.468e-04,3.565e-04,5.766e-04,8.070e-04,1.047e-03,1.296e-03,1.554e-03,1.820e-03,2.093e-03,2.372e-03,2.657e-03,2.947e-03,3.240e-03,3.535e-03,3.832e-03,4.128e-03,4.423e-03,4.715e-03,5.003e-03,5.284e-03,5.558e-03,5.821e-03,6.073e-03,6.312e-03,6.534e-03,6.739e-03,6.923e-03,7.085e-03,7.222e-03,7.330e-03,7.409e-03,7.454e-03,7.464e-03,7.435e-03,7.364e-03,7.248e-03,7.084e-03,6.869e-03,6.600e-03,6.273e-03,5.885e-03,5.432e-03,4.910e-03,4.317e-03,3.649e-03,2.901e-03,2.071e-03,1.154e-03,1.470e-04,9.543e-04,2.154e-03,3.455e-03,4.861e-03,6.378e-03,8.007e-03,9.754e-03,1.162e-02,1.361e-02,1.574e-02,1.799e-02,2.038e-02,2.291e-02,2.558e-02,2.840e-02,3.137e-02,3.449e-02,3.777e-02,4.121e-02,4.482e-02,4.859e-02,5.253e-02,5.664e-02,6.093e-02,6.539e-02,7.003e-02,7.486e-02,7.987e-02,8.507e-02,9.045e-02,9.602e-02,1.018e-01,1.077e-01,1.139e-01,1.202e-01,1.268e-01,1.335e-01,1.404e-01,1.476e-01,1.549e-01,1.624e-01,1.701e-01,1.780e-01,1.861e-01,1.944e-01,2.029e-01,2.116e-01,2.204e-01,2.295e-01,2.387e-01,2.481e-01,2.577e-01,2.675e-01,2.774e-01,2.875e-01,2.978e-01,3.082e-01,3.188e-01,3.295e-01,3.403e-01,3.514e-01,3.625e-01,3.738e-01,3.852e-01,3.967e-01,4.083e-01,4.200e-01,4.319e-01,4.438e-01,4.558e-01,4.679e-01,4.800e-01,4.923e-01,5.045e-01,5.169e-01,5.292e-01,5.416e-01,5.540e-01,5.664e-01,5.789e-01,5.913e-01,6.037e-01,6.161e-01,6.284e-01,6.407e-01,6.530e-01,6.652e-01,6.773e-01,6.894e-01,7.013e-01,7.132e-01,7.249e-01,7.366e-01,7.480e-01,7.594e-01,7.706e-01,7.817e-01,7.925e-01,8.032e-01,8.137e-01,8.240e-01,8.341e-01,8.440e-01,8.537e-01,8.631e-01,8.723e-01,8.812e-01,8.898e-01,8.982e-01,9.063e-01,9.141e-01,9.217e-01,9.289e-01,9.358e-01,9.423e-01,9.486e-01,9.545e-01,9.601e-01,9.653e-01,9.702e-01,9.748e-01,9.789e-01,9.827e-01,9.862e-01,9.892e-01,9.919e-01,9.942e-01,9.962e-01,9.977e-01,9.989e-01,9.997e-01,1.000e+00,1.000e+00,9.997e-01,9.989e-01,9.977e-01,9.962e-01,9.942e-01,9.919e-01,9.892e-01,9.862e-01,9.827e-01,9.789e-01,9.748e-01,9.702e-01,9.653e-01,9.601e-01,9.545e-01,9.486e-01,9.423e-01,9.358e-01,9.289e-01,9.217e-01,9.141e-01,9.063e-01,8.982e-01,8.898e-01,8.812e-01,8.723e-01,8.631e-01,8.537e-01,8.440e-01,8.341e-01,8.240e-01,8.137e-01,8.032e-01,7.925e-01,7.817e-01,7.706e-01,7.594e-01,7.480e-01,7.366e-01,7.249e-01,7.132e-01,7.013e-01,6.894e-01,6.773e-01,6.652e-01,6.530e-01,6.407e-01,6.284e-01,6.161e-01,6.037e-01,5.913e-01,5.789e-01,5.664e-01,5.540e-01,5.416e-01,5.292e-01,5.169e-01,5.045e-01,4.923e-01,4.800e-01,4.679e-01,4.558e-01,4.438e-01,4.319e-01,4.200e-01,4.083e-01,3.967e-01,3.852e-01,3.738e-01,3.625e-01,3.514e-01,3.403e-01,3.295e-01,3.188e-01,3.082e-01,2.978e-01,2.875e-01,2.774e-01,2.675e-01,2.577e-01,2.481e-01,2.387e-01,2.295e-01,2.204e-01,2.116e-01,2.029e-01,1.944e-01,1.861e-01,1.780e-01,1.701e-01,1.624e-01,1.549e-01,1.476e-01,1.404e-01,1.335e-01,1.268e-01,1.202e-01,1.139e-01,1.077e-01,1.018e-01,9.602e-02,9.045e-02,8.507e-02,7.987e-02,7.486e-02,7.003e-02,6.539e-02,6.093e-02,5.664e-02,5.253e-02,4.859e-02,4.482e-02,4.121e-02,3.777e-02,3.449e-02,3.137e-02,2.840e-02,2.558e-02,2.291e-02,2.038e-02,1.799e-02,1.574e-02,1.361e-02,1.162e-02,9.754e-03,8.007e-03,6.378e-03,4.861e-03,3.455e-03,2.154e-03,9.543e-04,1.470e-04,1.154e-03,2.071e-03,2.901e-03,3.649e-03,4.317e-03,4.910e-03,5.432e-03,5.885e-03,6.273e-03,6.600e-03,6.869e-03,7.084e-03,7.248e-03,7.364e-03,7.435e-03,7.464e-03,7.454e-03,7.409e-03,7.330e-03,7.222e-03,7.085e-03,6.923e-03,6.739e-03,6.534e-03,6.312e-03,6.073e-03,5.821e-03,5.558e-03,5.284e-03,5.003e-03,4.715e-03,4.423e-03,4.128e-03,3.832e-03,3.535e-03,3.240e-03,2.947e-03,2.657e-03,2.372e-03,2.093e-03,1.820e-03,1.554e-03,1.296e-03,1.047e-03,8.070e-04,5.766e-04,3.565e-04,1.468e-04,5.204e-05,2.398e-04,4.162e-04,5.812e-04,7.347e-04,8.766e-04,1.007e-03,1.126e-03,1.233e-03,1.329e-03,1.414e-03,1.488e-03,1.552e-03,1.605e-03,1.648e-03,1.681e-03,1.705e-03,1.719e-03,1.725e-03];
        
    else
        error('Unknown experiment');
    end
    
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


if v==0
    
    fa=fa(1:4:end);
    
    for i=1:length(fa)
        sart_tmp(i)=ssFLASH(fa(i),TR,T1,0,1);
    end
    sart=mean(sart_tmp);
    
    sart=sart*thk2/thk;
    return;
end

TotalSeg=100;

step=thk2/TotalSeg;

if v*TR>step
    
    sart_tmp =zeros(1,TotalSeg);
    
    for j=1:TotalSeg
        
        z=(j-1)*step;
        
        nRF=  ceil(z/v/TR);% number of RF that has been experienced by this spin
        
        sz=sz_init;
        for irep=1:nRF
            
            z0=z-(nRF-irep)*v*TR;
            fa0=interp1(pos,fa,z0);
            
            if irep==nRF
                sart_tmp(j) = sz*sin(fa0/180*pi);
            end
            sz=sz*cos(fa0/180*pi);
            
            
            sz=1+(sz-1)*exp(-TR/T1);
        end
        
    end
    
else
    step=v*TR;
    nRF = ceil(thk2/step);
    
    sart_tmp =zeros(1,nRF);
    sz=sz_init;
    
    for j=1:nRF
        
        z=(j-1)*step;
              
        fa0=interp1(pos,fa,z);
        
        sart_tmp(j) = sz*sin(fa0/180*pi);
        sz=sz*cos(fa0/180*pi);
        
        sz=1+(sz-1)*exp(-TR/T1);
        
    end
    
end
% below is same as
sart=mean(sart_tmp(:));

sart=sart*thk2/thk;


