function set_lineWidth(lw,h)

if ~exist('h','var')
    h = gcf;
end


ch = findobj(h,'Type','Axes','Tag','');

for i=1:length(ch)
    
    lines = findobj(ch(i),'Type','line','LineStyle','-');
    for j=1:length(lines)
      set(lines(j),'LineWidth',lw);
    end
  
    
    
end




