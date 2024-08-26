thk=0.2;
TR = 0.052/2;

VENC=4;
T1s=3.1;
T1f=3.1;
lambda=1;
FA=15;

flow_ana=zeros(3,nv);
flow_num=zeros(3,nv);
vmean=zeros(3,nv);
eflow_ana=zeros(3,nv);
FA=[15,25,35];
 exp_sel='fl_fq_retroZ_mb';

 
 flow_pattern='plug';
 pulse_profile='sinc';
 
  
tof_total(1)=1;

pv0=0.1;
tof_total(2)=(1-pv0+pv0*exp(1i*pi/8));

 [v_calc,sart_calc,v_calc_pc,sart_calc_pc,spvs]=Prep_PV_V4TOF2Sur_PC(TR,T1s,0.006,1.59,thk,FA(1),flow_pattern,pulse_profile,exp_sel);
 
 %%
 vall=0:0.1:1;
 pv=zeros(11,50);
 v=zeros(11,50);
 
 for j=1:50
 
     for i=1:length(vall)
 disp([i,j]);    
tof_total(2)=(1-pv0+pv0*exp(1i*pi*vall(i)/VENC))+(randn_white(1)+1i*randn_white(1))*0.01;
     
[pv(i,j),v(i,j)]=PV_V4PC_Hamilton(tof_total,v_calc,sart_calc,v_calc_pc,sart_calc_pc,spvs(1),VENC,lambda,flow_pattern,false);%% the TE and T2 are just made up; if the same for both flow and static, then doesn't matter
 
 end
 end
 figure;plot(std(pv.*v,[],2))
 %%
 for k=1:11
     
 flow=pv(k,:).*v(k,:);
 
 sd(k)=std(flow);
 
 
 end
 
 
 
 
 
 