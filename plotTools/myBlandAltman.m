function myBlandAltman(d1,d2,par)
% myBlandAltman(d1,d2,par)
% par can contains:
% sym:
% mode: relative or other
% xrng:
% sid:
% xl, yl;

if ~exist('par','var')
    par=[];
end

if ~isfield(par,'sym')
    
    sym={'o','>','<','s','*','d','v','x','p','^','h'};
end

if isfield(par,'mode') && strcmp(par.mode,'relative')
    y = 100*(d1-d2)./(d1+d2)*2;
else
    y = (d1-d2);    
end

sig=std(y);
mn=mean(y);

x=(d1+d2)/2;
if ~isfield(par,'xrng')
  xrng=[min(x)*0.5,max(x)*1.5];
else
    xrng=par.xrng;
end

if isfield(par,'sid')
 sid=par.sid;
else
 sid=ones(1,length(y));   
end

usid=unique(sid);
clr=lines(length(usid));
for i=1:length(usid)
    sel=sid==usid(i);  
    plot(x(sel),y(sel),sym{i},'Color',clr(i,:),'MarkerSize',8);
    hold on;
end

hold on;
plot(xrng,[1,1]*mn,'k-');

plot(xrng,[1,1]*(mn+1.96*sig),'k--');
plot(xrng,[1,1]*(mn-1.96*sig),'k--');


set(gca,'FontSize',12);

if isfield(par,'xl')
xlabel(par.xl);
end

if isfield(par,'yl')
ylabel(par.yl);
end

xlim(xrng);

%ymin=min([min(y),mn-1.96*sig]);
%ymax=max([max(y),mn+1.96*sig]);

%ymin=mn-1.96*sig;
%ymax=mn+1.96*sig;

fprintf('Coefficient of repeatability = %3.2e; %3.2e%%\n',1.96*sig,1.96*sig/mean(x)*100);

margin=0.5*sig;


set(gca,'TickLength',[0.03,0.03]);

%ylim([ymin-margin,ymax+margin]);


box on;
