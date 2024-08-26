
figure;
axis;
xlim([0,1]);
ylim([0,1]);
box on;
set(gca,'XTick',[-1,2],'YTick',[-1,2]);
t=tic;

instruction = {'tilt left','tilt right','tilt up','tilt down','shift left','shift right','shift out','shift in'}; 
start_time=30:60:600;

while 1
    current_time=toc(t);
    i=find(current_time-start_time>0,1,'last');
    if isempty(i)
        i=0;
    end
    disp(i);
    if i>length(instruction)
        break;
    end
  if i>0&&current_time-start_time(i)<3 && current_time-start_time(i)>0
    txt=instruction{i};
    clr='r';
  else
     txt=sprintf('%d',round(start_time(i+1)-current_time));
     clr='k';
  end
      
    hold off;
    h=text(0.4,0.5,txt,'FontSize',80,'Color',clr);
    pause(0.1);
    delete(h);
end

