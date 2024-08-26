function rdk_rename_ts

lMT = 2;
rMT = 1;

movefile(sprintf('ts_cluster%d.mat',lMT),'ts_hMT+L.mat');
movefile(sprintf('ts_cluster%d.fig',lMT),'ts_hMT+L.fig');

movefile(sprintf('ts_cluster%d.mat',rMT),'ts_hMT+R.mat');
movefile(sprintf('ts_cluster%d.fig',rMT),'ts_hMT+R.fig');



