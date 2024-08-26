function set_fontSize(fs,h)

if ~exist('h','var')
    h = gcf;
end


ch = findobj(h,'Type','axes','Tag','');

for i=1:length(ch)
    
    set(ch(i),'FontSize',fs);
    
    text_ch = findobj(ch(i),'Type','text');
    
    for j=1:length(text_ch)
        set(text_ch(j),'FontSize',fs);
    end
    
    
set(get(ch(i),'xlabel'),'FontSize',fs);

set(get(ch(i),'ylabel'),'FontSize',fs);
set(get(ch(i),'title'),'FontSize',fs);
    
end




