s=[1131,1481.6,1808.8,2046.3,2140.4,2143.2,2042.4,1840.0,1592.7,1270.1];
B1=90:20:270;
%% measured Ref ampltidue = 190.8 V.

figure;plot(B1,s,'o-');

set(gca,'FontSize',14);

xlabel('Reference Amplitude (V)');
ylabel('Image Intensity (arb. units)');

xlim([60,280]);