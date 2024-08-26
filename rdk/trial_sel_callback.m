function trial_sel_callback(a,nTR_trial)


cp = get(gca,'CurrentPoint');

it = floor(cp(1,1)/nTR_trial);
ch_r = findobj(gca,'Type','Rectangle');

find_rect = false;
t_exc = [];
for i=1:length(ch_r)
    rect_pos = get(ch_r(i),'Position');
    if rect_pos(1)/nTR_trial == it
        delete(ch_r(i));
        find_rect = true;
    else
        t_exc(end+1) = rect_pos(1)/nTR_trial+1;
    end
end

if ~find_rect
 yl=ylim;
 r=rectangle('Position',[it*nTR_trial,yl(1),nTR_trial,diff(yl)*0.2]);
 set(r,'FaceColor',[0.5,0.5,0.5]);
 t_exc(end+1) = it+1;
end
   

% plot the mean and sem.

l=findobj(gca,'Type','line');
for i=1:length(l)
    x=get(l(i),'XData');
    if x(2)-x(1) == 1;
     y = get(l(i),'YData');
     break;
    end
end

fprintf('[%s]\n',num2str(sort(t_exc)));

d = reshape(y,nTR_trial,length(y)/nTR_trial);
d(:,t_exc) = [];
 mn = mean(d,2);
 sem = std(d,0,2)/sqrt(size(d,2));

 fprintf('total variance: %4.3f\n',sum(sem.^2)*10000);
cf = gcf; 
figure(a);  errorbar(mn,sem,'r');
figure(cf);
