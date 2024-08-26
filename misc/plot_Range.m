function [xrange,yrange]=plot_Range(handle,type)
% [xrange,yrange]=plot_Range(handle,type)
% handle: handle of the parent of plotted data
% type: data type.
ch= findobj(handle,'Type',type);
xdata = [];
ydata = [];
for i=1:length(ch)
    d_temp=get(ch(i),'XData');
    xdata = [xdata;d_temp(:)];
    d_temp = get(ch(i),'YData');
    ydata = [ydata;d_temp(:)];
end

xrange = zeros(1,2);
yrange = zeros(1,2);
xrange(1) = min(xdata(:));
xrange(2) = max(xdata(:));

yrange(1) = min(ydata(:));
yrange(2) = max(ydata(:));