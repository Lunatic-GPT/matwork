function savetiff(figname,handle)
% save to tiff using print() function.
% matlab figure resolution 90dpi
% screen resolution 150 dpi. Therefore, the fonts are smaller than they
% appear on screen.
% savetiff(figname,handle,width)
% width is in inch.
if ~exist('handle','var') || isempty(handle)
    handle = gcf;
end


figure(handle);

set(gca,'XTickMode','manual','YTickMode','manual');
%set(gca,'XTickLabelMode','manual','YTickLabelMode','manual');
set(gca,'XLimMode','manual','YLimMode','manual');
set(gcf,'Units','inches');
pos = get(gcf,'Position');

set(gcf,'PaperUnits','inches');
set(handle,'PaperPositionMode','manual');

pos(1) = 0;
pos(2) = 0;
set(handle,'PaperPosition',pos);

prefix = strrm(figname,'.tif');
saveas(handle,[prefix,'.fig']);
pause(0.1);
print(handle,'-r300','-dtiff',[prefix,'.tif']);

%cropimage([prefix,'.tif']);



