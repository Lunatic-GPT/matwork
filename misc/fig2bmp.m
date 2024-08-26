function fig2bmp(pattern)

str = dir(pattern);
for i=1:length(str)
  name = strtok(str(i).name,'.');    
  open(str(i).name);
  %set(gcf,'Units','inch','Position',[2   3   2.5   1.8]);
  set(gca,'Box','off');
   %ch = get(gcf,'Children');
   %set(ch(2),'Position',[0.5703    0.18    0.33    0.70]);
   %set(ch(4),'Position',[0.1300    0.18    0.33    0.70]);
  %ch = get(gcf,'Children');
  %set(ch(2),'Box','off');
  %set(ch(2),'TickLength',[0.02,0]);
  
  saveas(gcf,[name,'.png']);
  close;
  if ~exist('bmp_fig','dir')
     mkdir('bmp_fig');
  end
  movefile('*.png','bmp_fig');
end
