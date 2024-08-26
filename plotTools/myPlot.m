function h=myPlot(x,y,symb,par)
% myPlot(x,y,symb,par)
% x is n*1 or n*m;
% y is n*1 or n*m;

if size(x,1)==1
    x=x(:);
end
if size(y,1)==1
    y=y(:);
end

if ~exist('par','var') || ~isfield(par,'color') 
    Color = lines(size(y,2));
else
    Color=repmat(par.color,[ceil(size(y,2)/size(par.color,1)),1]);
end

if size(x,2)==1 && size(y,2)>1
    x=repmat(x,[1,size(y,2)]);
end

symb=str2cell(symb);
symb=repmat(symb,[1,size(y,2)]);

for i=1:size(y,2)

    h(i)= plot(x(:,i),y(:,i),symb{i},'LineWidth',1);
    hold on;
    set(gca,'FontSize',12);
    par.color=Color(i,:);
    if exist('par','var')
      set_plot(h(i), par);
    end

end


