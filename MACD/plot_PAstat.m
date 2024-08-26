function plot_PAstat(mbac,clr,sym)

mbac=str2cell(mbac);

for i=1:length(mbac)
load(mbac{i});
vall(:,i)=v';
radall(:,i)=rad';
end

inc=vall>0&radall<0.2;
nv=size(vall,1);
par.xlabel='Vessel';
par.ylabel='V (cm/s)';
par.color=clr;
subplot(1,2,1);
myPlot(1:size(vall,1),vall,sym,par);
hold on;
myPlot(1:size(vall,1),mean(vall(inc))*ones(1,nv),['-',clr],par);

subplot(1,2,2);

par.ylabel='D (mm)';
myPlot(1:size(vall,1),radall*2,sym,par);
hold on;
myPlot(1:size(radall,1),mean(radall(inc))*2*ones(1,nv),['-',clr],par);




