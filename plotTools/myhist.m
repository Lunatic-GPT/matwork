function myhist(d,xbin,xl,yl)

hist(d,xbin);
   set(gca,'FontSize',12);
   xlim(xbin([1,end-1]));
   
   xlabel(xl);
   ylabel(yl);
   
   c=get(gca,'Children');
   
   set(c,'FaceColor',0.5*[1,1,1]);
%   set(c,'FaceAlpha',0);
   set(gca,'TickLength',[0.02,0.02],'TickDir','out');
   