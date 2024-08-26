function axis_scale(height)

ch = findobj(gcf,'Type','axes','Tag','');

for i=1:length(ch)
    pos = get(ch(i),'Position');
   
    pos(4) = height;
    set(ch(i),'Position',pos);
end
