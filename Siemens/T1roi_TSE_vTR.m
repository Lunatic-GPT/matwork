function T1roi_TSE_vTR(dname,roi)
% T1roi_TSE_vTR(dname,roi)
d=ri(dname,1);

extp(dname);
TRmin=readsPar([dname,'.pro'],'adFree[5]');

TRmax=readsPar([dname,'.pro'],'adFree[4]');

nTR = readsPar([dname,'.pro'],'alFree[5]');

logsp=readsPar([dname,'.pro'],'alFree[4]');

if logsp
    TR=logspace(log10(TRmin),log10(TRmax),nTR);
else
    TR = linspace(TRmin,TRmax,nTR);
end


roi=ri(roi);
y=mean_roi(d,roi);

figure;plot(TR,y,'o');

fittool(@exp_decay2,[-max(y),mean(TR(1:5)),max(y)]);

