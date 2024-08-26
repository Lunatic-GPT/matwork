a=load('TR_0.03_T1_2.6_1.2_T2_0.0285_0.0235_TE_0.0157_thk_0.2_FA_36_tV_dvdt-2to2.mat');
b=load('TR_0.03_T1_2.6_1.2_T2_0.0285_0.0235_TE_0.0157_thk_0.2_FA_36_tV_dvdt-1to1.mat');


v_pc=a.v_pc;
dvdt=[a.dvdt,b.dvdt];
s_pc=[a.s_pc;b.s_pc];

[dvdt,ind]=sort(dvdt);
s_pc=s_pc(ind,:);
s_static=a.s_static;

save('TR_0.03_T1_2.6_1.2_T2_0.0285_0.0235_TE_0.0157_thk_0.2_FA_36_tV.mat','v_pc','dvdt','s_pc','s_static');


%%


