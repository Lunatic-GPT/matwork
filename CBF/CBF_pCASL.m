function res=CBF_pCASL(sl,sc,T1a,eff,tlab,tran,T1ts,T1t,pc)
% s1: label signal
% sc: control signal
% T1a: artery T1 in s
% eff: labelling efficiency
% tlab: labeling time
% tran: transit time
% T1ts: tissue T1 including CBF contribution
% T1t: tissue T1 without CBF contribution
% pc: tissue-blood partition coefficient

if ~exist('pc','var')
    pc=0.9; %unit ml/g
end
eff=eff*exp(-tran/T1a);

eta=(1-exp(-(tlab-tran)/T1ts));
ds=sc-sl;

res=pc/T1t*ds./(2*eff*sc*eta-ds);
