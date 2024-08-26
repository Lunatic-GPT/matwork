function plot_PC_CardiacCycle(dcm,mask,mask_bg,venc)
d=ri(dcm);
m=ri(mask);
d2=mean_roi(d,m);
mb=ri(mask_bg);
mb=mb&~m;
ph0=mean_roi(mean(d,4),mb);
v=(d2-ph0)*venc/2048;
x=(0:size(d,4)-1)/size(d,4);

figure;plot(x,v,'o-');
set(gca,'FontSize',12);
xlabel('Caridac cycle');
ylabel('Apparent velocity (cm/s)');