function myboxplot(xtklb,ylbl,ytk,varargin)
% myboxplot(xtklb,ylbl,ytk,varargin)
x=[];
group=[];

if nargin>4
    for i=4:nargin
        x=cat(1,x,varargin{i-3}(:));
        group=cat(1,group,i*ones(length(varargin{i-3}),1));
    end
elseif iscell(varargin{1})  % cell array
    for i=1:length(varargin{1})
        x=cat(1,x,varargin{1}{i});
        group=cat(1,group,i*ones(length(varargin{1}{i}),1));
    end
else  % 2D array;
     for i=1:size(varargin{1},2)
        x=cat(1,x,varargin{1}(:,i));
        group=cat(1,group,i*ones(size(varargin{1},1),1));
    end
end

boxplot(x,group);
set(gca,'FontSize',12);
set(gca,'XTickLabel',xtklb);
ylabel(ylbl);

if ~isempty(ytk)
set(gca,'YTick',ytk);
ylim(ytk([1,end]));
end
hm=findobj(gca,'tag','Median');

ho=findobj(gca,'tag','Outliers');

hLine=findobj(gca,'Type','Line');

    set(hLine,'Color','k','LineWidth',1);
    set(hm,'Color','b');
    
    set(gca,'TickLength',[0.02,0.02]);
    set(ho,'MarkerEdgeColor','b','MarkerSize',10,'Marker','.');
    
box on;