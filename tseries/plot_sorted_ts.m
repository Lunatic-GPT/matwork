function [mn,sem]=plot_sorted_ts(ts,nTR_trl)

%plot_sorted_ts(ts,nTR_trl)
  d = ts;
  [mn,sem] = ts_sort2ave1d(ts,nTR_trl);
  h=figure;
  %errorbar(mn,sem,'r');
  plot(mn,'r');
  figure;
    plot(0:length(d)-1,d);
    xlim([0,length(d)]);
 hold on;
  ind = nTR_trl:nTR_trl:length(d(:));
 plot(ind-1,d(ind),'or');
  ind = 0:nTR_trl:length(d)-nTR_trl;
 plot(ind,d(ind+1),'*k');
 set(gca,'XGrid','on');
 set(gca,'XTick',0:nTR_trl:length(d(:)));
 
 set(gca,'ButtonDownFcn',sprintf('trial_sel_callback(%d,%d)',h,nTR_trl));
 

