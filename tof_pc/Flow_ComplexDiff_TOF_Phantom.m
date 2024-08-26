function [flow_num,pv,v,water_fraction]=Flow_ComplexDiff_TOF_Phantom(flist_tof,flist_pc,FA,TR,T1,thk,VENC,m_bg,m_tube,water_fraction,vox_size,neg_phase,exp_sel,twoPool,flow_pattern,pulse_profile)
       
% flist_tof: the list of tof scans.  The first one should have no flow to
% correct for the signal void of tube.
% flist_pc: same as convention and length as flist_tof.  PC should be in
% units of degree
% water_fraction: the water fraction in m_tube
% 2/24/2017: 1. add phase and tof within the roi first before calculate pv and v maps
%            2. tof calculated wrt surrounding static tissue instead of wrt
%            the same voxel without flow.
% venc: cm/s
% vox_size: cm2
% flow: cm3/s
% vmean: cm/s
% thk: cm
% err: error in the magnitude image

%  flow_pattern='plug';
%  pulse_profile='sinc';
  
m_pvf1=ri(m_bg,'','','d');

%m=repmat2(m_tmp,10);
m=ri(m_tube,'','','d');

if length(FA)==1
    [v_calc,sart_calc,v_calc_pc,sart_calc_pc,s_static]=Prep_PV_V4TOF2Sur_PC(TR,T1,0.006,1.59,thk,FA,flow_pattern,pulse_profile,exp_sel);
end

s_static=s_static(1);

% tof_noflow=ri(flist_tof{1});
% tof_noflow=tof_noflow(:,:,:,1);
% 
% %tof_noflow=tof_noflow(1:10:end,1:10:end);
% 
% tof_noflow_full=mean(tof_noflow(m_pvf1>0));
%   

  nv=sum(m(:)>0);
%   water_fraction = sum(tof_noflow(m>0))/tof_noflow_full/nv;
%   
nf=length(flist_tof);
ok=zeros(1,nf);
%flow_nocorr=zeros(1,nf);
pv=zeros(1,nf);
v=zeros(1,nf);

disp(['tof    ','pv    ','v    ','pv*v    ']);

for j=1:length(flist_tof)
    try
       pc=ri(flist_pc{j},[],[],'d');
    catch
        pc=ri(flist_pc{j});
    end
   if neg_phase
       pc=-pc;
   end
   
  tof_flow=ri(flist_tof{j});
  
  %pc=pc(1:10:end,1:10:end);
  %tof_flow=tof_flow(1:10:end,1:10:end,:,:);
  
  tof_flow_nofe=tof_flow(:,:,1,1);
  tof_flow_fe=tof_flow(:,:,1,2);
  
  tof_full_nofe=mean(tof_flow_nofe(m_pvf1>0));
  %tof_full_fe=mean(tof_flow_fe(m_pvf1>0));  % should we normalize them separately?
  tof_full_fe=tof_full_nofe;  % should we normalize them separately?
  
  
  %water_faction=0.9906;
  
  pc_bg=mean(pc(m_pvf1>0));
  
  M_total(1)=sum(tof_flow_nofe(m>0))/tof_full_nofe/nv/water_fraction;
  M_total(2)=sum(tof_flow_fe(m>0).*exp(1i*(pc(m>0)-pc_bg)/180*pi))/tof_full_fe/nv/water_fraction;
  
  z(j,:)=[M_total(2)-M_total(1), M_total(1)];
  
  %pc_total=sum((pc(m>0)-pc_bg).*tof_noflow(m>0))/sum(tof_noflow(m>0));
  
  
  lambda=1;
  verbose=true;
  if length(FA)>1
    [v_calc,sart_calc,v_calc_pc,sart_calc_pc,s_static]=Prep_PV_V4TOF2Sur_PC(TR,T1,0.006,1.59,thk,FA(j),flow_pattern,pulse_profile,exp_sel);%% the TE and T2 are just made up; if the same for both flow and static, then doesn't matter
  end
  
  s0=sart_calc(1);
  sart_calc=sart_calc/s0;
  sart_calc_pc=sart_calc_pc/s0;
  s_static=s_static/s0;
  
  
  [pv(j),v(j),ok(j)]=PV_V4ComplexDiff_TOF(z(j,:), v_calc,sart_calc,v_calc_pc,sart_calc_pc,s_static,VENC,lambda,flow_pattern,twoPool);
    
  
  
  %flow_nocorr(j)=sum(pc(m>0).*tof_noflow(m>0))/tof_noflow_full*VENC/180*vox_size;
  
  
  %%
end


%disp(tof_total);
%disp([pv,v,pv*v]);
  flow_num=pv.*v*vox_size*nv*water_fraction;  %cm3/s  
 

