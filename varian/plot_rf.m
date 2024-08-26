function plot_rf(fname)
% plot_rf(fname)

a=textread(fname,'','commentstyle','shell');

%rf=a(:,2).*exp(1i*a(:,1)/180*pi);

freq = diff(a(:,1)/360)*size(a,1);
figure;subplot(1,2,1);plot(-freq);
title('Frequency (unit 1/pw)');
xlabel('Time');

subplot(1,2,2);plot(abs(a(:,2)));
title('Amplitude');
xlabel('Time');

