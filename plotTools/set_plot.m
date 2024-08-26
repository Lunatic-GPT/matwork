function set_plot(h,par)

if isfield(par,'marker')
    set(h,'Marker',par.marker);
end


if isfield(par,'markersize')
    set(h,'MarkerSize',par.markersize);
end

if isfield(par,'xlim')
    xlim(par.xlim);
end

if isfield(par,'fontsize')
    set(gca,'FontSize',par.fontsize);
end

if isfield(par,'linewidth')
    set(h,'LineWidth',par.linewidth);
end

if isfield(par,'linestyle')
    set(h,'LineStyle',par.linestyle);
end

if isfield(par,'ylim') && ~isempty(par.ylim)
    ylim(par.ylim);
end

if isfield(par,'ylabel')
    ylabel(par.ylabel);
end

if isfield(par,'xlabel')
    xlabel(par.xlabel);
end

if isfield(par,'color')
    set(h,'Color',par.color);
end

if isfield(par,'title') 
    title(par.title);
end

if isfield(par,'legend') 
    legend(par.legend);
end

set(gca,'TickLength',[0.02,0.02]);

if isfield(par,'xtick')
    set(gca,'XTick',par.xtick);
end

if isfield(par,'xticklabel')
    set(gca,'XTickLabel',par.xticklabel);
end

if isfield(par,'ytick')
    set(gca,'YTick',par.ytick);
end


 box on;
