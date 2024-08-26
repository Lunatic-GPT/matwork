function res=CBF_pCASL_Artery(ds,tlab,tpostlab,pc,T1a,eff)
% res=CBF_pCASL_Artery(ds[,tlab,tpostlab,pc,T1a,eff])
% ds: normalized signal change (sc-sl)/sc
% tlab: labeling time
% tpostlab: postlab delay
% pc: tissue-blood partition coefficient
% T1a: artery T1 in s
% eff: labelling efficiency
% result in ml/100g/min
if ~exist('pc','var')
    pc=0.9; %unit ml/g
end

if ~exist('eff','var')|| isempty(eff)
    eff=0.68;
end

if ~exist('T1a','var')|| isempty(T1a)
    T1a=1/0.67; %unit s
end

if ~exist('tpostlab','var')|| isempty(tpostlab)
    tpostlab=1; %adFree[2]/1000000
end

if ~exist('tlab','var')|| isempty(tlab)
    tlab=18.5*82/1000; % adFree[3]*18.5/1000
end


res=pc*ds/T1a/2/eff/exp(-tpostlab/T1a)/(1-exp(-tlab/T1a));

res=res*6000;%% convert to ml/100g/min

