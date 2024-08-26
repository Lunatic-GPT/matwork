function plot_motion_mat(fmat,TR,voxsize)

t=load(fmat);


time=TR*(0:size(t.tran,1)-1);
figure;
clr='rbg';
    subplot(1,2,1);
    for i=1:3

        if i==2
            y=-t.tran(:,i)*voxsize;
        else
             y=t.tran(:,i)*voxsize;
        end
        plot(time,y,'.-','Color',clr(i));
        hold on;
        if i==3
            xlabel('Time (s)');
        end
        
    end
    ylabel('Translation (mm)');
    xlim([0,max(time)+1]);
    title('Translation');
  %  legend(lbl(1:3));
    set(gca,'FontSize',12);
    subplot(1,2,2);
    clr='brg';
    for i=1:3
        
        if i==2
            y=-t.angle(:,i)*voxsize;
        else
             y=t.angle(:,i)*voxsize;
        end

        plot(time,y,'-','Color',clr(i));
        hold on;

        if i==3
            xlabel('Time (s)');
        end
        
    end
    
    ylabel('Rotation (deg)');
    xlim([0,max(time)+1]);
    title('Rotation');
%    legend(lbl(1:3));
    
    set(gca,'FontSize',12);

