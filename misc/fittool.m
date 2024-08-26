function beta=fittool(fun,b0,fit_xrange)
%fittool(fun,b0[,fit_xrange])
if strcmp(get(gco,'Type'),'line')
   h = gco;
else
  h = findobj(gca,'Type','line');
  
  for i=1:length(h)
      clr=get(h(i),'Color');
      
      if any(clr~=[1,0,0])
          h=h(i);  % find the first curve that is not red (i.e., not fitted result)
          break;
      end
  end
end

x=get(h,'XData');
y=get(h,'YData');

if ~exist('fit_xrange','var')
  fit_xrange=min_max(x);
end
xl = fit_xrange(1);
xh = fit_xrange(2);

ind = find(x>=xl & x<=xh);
%options=optimset('FunValCheck','off');
options=optimset('MaxFunEvals',150);
beta=lsqcurvefit(fun,b0,x(ind),y(ind),[],[],options);
%[beta,resnorm,resi,exitflag,ouput]=lsqcurvefit(fun,b0,x(ind),y(ind));

for i=1:length(beta)
    fprintf('Parameter %d = %f\n',i,beta(i));
end

hold on;
%xx = linspace(xl,xh,1000);
xx=sort(x);
plot(xx,fun(beta,xx),'r-');
mxx = median(xx);
a=help(func2str(fun));
format = textscan(a,'%s');
text(mxx,fun(beta,mxx),sprintf(format{1}{end},beta),...
     'HorizontalAlignment','left','FontName','Arial');




