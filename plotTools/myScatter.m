function [b,h]=myScatter(x,y,draw_line,par,par_fit)
% myScatter(x,y,draw_line[,par,par_fit])
% draw_line: 0-no line; 
%            1-draw line of indentity
%            2-draw regression line (default);
%            3 - draw regression line with intercept = 0

    sym='k.';
x=double(x);
y=double(y);
if ~exist('draw_line','var')
    draw_line=2;
end

if ~exist('par_fit','var')
    par_fit=[];
end

if ~exist('par','var')
    par=[];
end


h=plot(x,y,sym,'MarkerSize',6,'LineWidth',1);
set(gca,'FontSize',12);

if draw_line>0
disp_corr(x,y);

disp_corr(x,y,'Spearman');
end
b=0;
clr=lines(3);
if draw_line==1
    hold on;
    mx=plot_range([x(:);y(:)],0.1);
    
%     if mx(1)<0 && mx(2)<0
%         plot([0,mx(1)],[0,mx(1)],'k-');
%     elseif mx(1)<0 && mx(2)>0
%         
   h(2)=     plot(mx,mx,'-','LineWidth',1,'Color',clr(2,:));
%     else
%         plot([0,mx(2)],[0,mx(2)],'k-');
%     end
    
elseif draw_line==2
    hold on;
    
    xx=[x(:),ones(length(x),1)];
    
    b=xx\y(:);
    fprintf('b (y-intercept) = %f; k=%f; x-intercept = %f\n', b(2),b(1),-b(2)/b(1));
    xx2=[plot_range(x,0.1);1,1]';
    
  h(2)=  plot(xx2(:,1),xx2*b,'-','LineWidth',1,'Color',clr(2,:));
    
    xlim(plot_range(x,0.1));
    ylim(plot_range(y,0.1));
    
elseif draw_line==3
    
     hold on;
    
    xx=x(:);
    
    b=xx\y(:);
    xx2=[plot_range(x,0.1)]';
    
   h(2)= plot(xx2(:,1),xx2*b,'-','LineWidth',1,'Color',clr(2,:));
    
    xlim(plot_range(x,0.1));
    ylim(plot_range(y,0.1));
    
    fprintf('k=%f\n', b(1));
end

set(gca,'TickLength',[0.02,0.02]);

if draw_line==0
for i=1:length(h)
    set_plot(h(i),par);
end
else
set_plot(h(1),par);


set_plot(h(2),par_fit);
end



