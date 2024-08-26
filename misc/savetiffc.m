function savetiffc(figname,hand)
% save to tiff using print() function.
% matlab figure resolution 90dpi
% screen resolution 150 dpi. Therefore, the fonts are smaller than they
% appear on screen.
% savetiff(figname,hand,width)
% width is in inch.
if ~exist('hand','var') || isempty(hand)
    hand = gcf;
end

figure(hand);

set(gca,'XTickMode','manual','YTickMode','manual');
%set(gca,'XTickLabelMode','manual','YTickLabelMode','manual');
set(gca,'XLimMode','manual','YLimMode','manual');
set(gcf,'Units','inches');
pos = get(gcf,'Position');

set(gcf,'PaperUnits','inches');
set(hand,'PaperPositionMode','manual');

pos(1) = 0;
pos(2) = 0;
set(hand,'PaperPosition',pos);

prefix = strrm(figname,'.tif');
saveas(hand,[prefix,'.fig']);
print('-r300','-dtiff',[prefix,'.tif']);

cropimage([prefix,'.tif']);
delete([prefix,'.tif'])



