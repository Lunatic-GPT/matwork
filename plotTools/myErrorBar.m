function h=myErrorBar(x,y,e,par)
%  h=myErrorBar(x,y,e,par)
markerSize=6;
if ~exist('par','var')
    par=struct;
end
if ~isfield(par,'color') 
    color = lines(size(y,2));
else
    color=par.color;
end

if ~isfield(par,'marker')
    par.marker='s';
end

if isfield(par,'marker')&& length(par.marker)==1
    marker = repmat(par.marker,[1,size(y,2)]);
else
    marker=par.marker;
end


if ~isfield(par,'fontsize')
    par.fontsize=12;
end

if ~isfield(par,'linestyle')
    linestyle='-';
else
    linestyle=par.linestyle;
end

%set(gca,'FontSize',12);

if (size(x,1)==1||size(x,2)==1) && (size(y,1)>1&&size(y,2)>1) %x 1*n or n*1, y: n*m;
    
    for i=1:size(y,2)
        par.color=color(i,:);
        par.marker=marker(i);
        
        h(i)= errorbar(x(:),y(:,i),e(:,i),'LineWidth',1,'MarkerSize',markerSize,'LineStyle',linestyle);
        set_plot(h(i),par);
        hold on;
    end
elseif ~any(size(x)~=size(y)) &&(size(y,1)>1&&size(y,2)>1) % x:n*m, y:n*m

    for i=1:size(y,2)
        par.color=color(i,:);
        par.marker=marker(i);
        
        h(i)= errorbar(x(:,i),y(:,i),e(:,i),'LineWidth',1,'MarkerSize',markerSize,'LineStyle',linestyle);
        set_plot(h(i),par);
        hold on;
    end 
else % x,y both vectors
    h=errorbar(x,y,e,'LineWidth',1,'MarkerSize',markerSize,'LineStyle',linestyle);
     set_plot(h,par);
end



box on;
