function calc_CBF_pCASL(dname)

extract_protocol(dname);

a=ri(dname);

T1t=1.4;
T1a=1.8;
eff=0.9;

res=CBF_pCASL(sl,sc,T1a,eff,tlab,tran,T1ts,T1t);

save(['CBF_',dname],'res');