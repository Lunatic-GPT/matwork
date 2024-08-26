function plot_stimpatch(interval)

% interval is a n*2 matrix

yl=ylim;
for i=1:size(interval,1)
    
    hold on;
    h=patch([interval(i,1)*ones(1,2),interval(i,2)*ones(1,2)],[yl,yl(2),yl(1)],0.6*ones(1,3));
    set(h,'FaceAlpha',0.5)
end
