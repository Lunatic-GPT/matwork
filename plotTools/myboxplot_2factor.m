function myboxplot_2factor(par,varargin)
% each varargin{i} should contain two or more vectors that will be placed
% close together for comparison.
x=[];
group1=[];
group2=[];

if nargin>2
ng1=length(varargin{1});
ng2=nargin-1;

    for i=1:ng2
        for j=1:ng1
          x=cat(1,x,varargin{i}{j}(:));
          group1=cat(1,group1,j*ones(length(varargin{i}{j}),1));
          group2=cat(1,group2,i*ones(length(varargin{i}{j}),1));
        end
    end
else
ng1=size(varargin{1},1);
ng2=size(varargin{1},2);

    for i=1:ng2
        for j=1:ng1
          x=cat(1,x,varargin{1}{j,i}(:));
          group1=cat(1,group1,j*ones(length(varargin{1}{j,i}),1));
          group2=cat(1,group2,i*ones(length(varargin{1}{j,i}),1));
        end
    end

end
boxplot(x,[group2,group1],'factorseparator',1,'factorgap',10);

set(gca,'FontSize',12);


hm=findobj(gca,'tag','Median');

ho=findobj(gca,'tag','Outliers');

hLine=findobj(gca,'Type','Line');

    set(hLine,'Color','k','LineWidth',1);
    set(hm,'Color','b');
    
    set(gca,'TickLength',[0.02,0.02]);
    set(ho,'MarkerEdgeColor','b','MarkerSize',10,'Marker','.');
    
box on;

set_plot(gca,par);

