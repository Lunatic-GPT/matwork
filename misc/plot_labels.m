function plot_labels(title_str,xlab,ylab,leg,fs)
%plot_labels(title_str,xlab,ylab,leg,fs)


 title(title_str,'FontSize',fs);
 xlabel(xlab,'FontSize',fs);
 ylabel(ylab,'FontSize',fs);
 if isempty(leg)
     legend off;
 else
 legend(leg,'FontSize',fs-2);
 end
 set(gca,'FontSize',fs);
 
 
 