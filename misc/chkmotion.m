function chkmotion(pc)
% chkmotion(pc[=0])
% pc: 1 for physiological corrected, 0 otherwise
if ~exist('pc','var')
    pc = 0;
end

a = dir('motion_s*.1D');
for i=length(a):-1:1
    if ~pc
      d=load(sprintf('motion_s%d.1D',i));
    else
        d=load(sprintf('pcmotion_s%d.1D',i));
    end
    figure;
    set(gcf','Units','inch','Position',[3,3,8,10]);
    tit = {'A-P[mm]','R-L[mm]','I-S[mm]','Yaw[o]','Pitch[o]','Roll[o]'};
    for j=1:size(d,2)
        subplot(6,1,j);
        plot(0:size(d,1)-1,d(:,j),'-');
        legend(tit{j},'Location','EastOutside');
        if j==1
            title(sprintf('Scan %d',i));
        end
    end
  %  saveas(gcf,sprintf('motion_scan%d.tif',i));
end