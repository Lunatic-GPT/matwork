x=linspace(1,100);

figure;
%% snr vs scan time
subplot(1,2,1);
plot(x,1./sqrt(x),'b','LineWidth',2);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
xlim([0,110]);

ylim([0,1.1]);

%% snr vs voxel size
subplot(1,2,2);
plot(x,1./(x),'b','LineWidth',2);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
xlim([0,100]);

ylim([0,1.2]);
