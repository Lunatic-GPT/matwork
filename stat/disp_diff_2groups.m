function disp_diff_2groups(x,y,type,par)
% type can be:
% pt: paired T
% ut: un-paired T
% sr: signed-rank
% rs: rank sum
 par.linestype='-';
 par.xtick=[1,2];
par.xlim=[0.8,2.2];

if strcmp(type,'pt')
    p=ttest_1vec(x-y);
    myPlot([1,2],[x(:),y(:)]','o-',par);
elseif strcmp(type,'ut')
    p=ttest_2vec(x,y);
    myPlot([1,2],[x(:),y(:)]','o',par);
elseif strcmp(type,'sr')
    p=signrank(x-y);    
    myPlot([1,2],[x(:),y(:)]','o-',par);
elseif strcmp(type,'rs')
    p=ranksum(x,y);
    myPlot([1,2],[x(:),y(:)]','o',par);
else
    error('unknown type');
end
xl=xlim;
yl=ylim;
xpos=xl(1)+0.2*(xl(2)-xl(1));
ypos=yl(2)-0.2*(yl(2)-yl(1));

    text(xpos,ypos,sprintf('p = %3.2g',p),'FontSize',12);
fprintf('group difference: %g vs %g (%g%%); \n',mean(x),mean(y),mean(y./x-1)*100);
fprintf('p value = %g\n',p);



