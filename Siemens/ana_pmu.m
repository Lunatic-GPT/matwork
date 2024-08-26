a=load('Physio_fl_fq_retro_noheader.txt');

[actualData,ushLine,ushPartition,ushSlice,ushSet,ulPMUTimeStamp,freePara] = readMeasDat...
    ('meas_MID14_fl_fq_retroZ_sb_Sinc_PAT2_105V_pilot_FID14444.dat',inf,0,true);


figure;plot(a(5:end));

%%
mdhstart=57082387;
mdhstop=57189337;

atmp=a;
atmp(a>4096)=[];
figure;plot((0:length(atmp)-1)*0.02,atmp);
%%

n=length(a)-4-length(find(a>4096));

dt=mdhstop-mdhstart;

a2=a(5:end-1);
a2=a2(5093:9334);
a2=a2(28:end);

a2(a2>4095)=[];%2048;
t0=mdhstart;
t=linspace(mdhstart,mdhstop,length(a2));

figure;%plot((t-t0)/1000,a2);
plot((0:length(a2)-1)*0.02,a2);  % the log file sampling rate is 20 ms

hold on;
t0b=min(ulPMUTimeStamp);  % the PMUTimeStamp in MDH is in units of 2.49 ms

amdh=freePara(1:32:end,4);
amdh=amdh(82:end);

ts=ulPMUTimeStamp(1:32:end);
ts=ts(82:end);
ts=ts-ts(1);
%plot((0:length(amdh)-1)*0.026,amdh,'r');

plot(ts*0.00249,amdh,'r');

