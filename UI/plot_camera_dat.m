
function plot_camera_dat(dname)

fname=name4pat(fullfile(dname,'*.3dcamera.dat'),1);
fid = fopen(fname);

C = textscan(fid, '%s',23);
d = textscan(fid, '%f,%f,%f,%f,%f,%f,%f',Inf);

fclose(fid);

      %%
figure;
subplot(1,2,1);
plot(d{1}/100000,d{2},'r-');
hold on;
plot(d{1}/100000,d{3},'g-');
plot(d{1}/100000,d{4},'b-');
xlabel('Time (s)');
ylabel('Translation (mm)');
title('Translation');
legend('X','Y','Z');

subplot(1,2,2);
plot(d{1}/100000,d{5}*180/pi,'r-');
hold on;
plot(d{1}/100000,d{6}*180/pi,'g-');
plot(d{1}/100000,d{7}*180/pi,'b-');
xlabel('Time (s)');
ylabel('Rotation (deg)');
title('Rotation');
legend('X','Y','Z');









