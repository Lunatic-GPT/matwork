function plot_motion_afni(dfile,separate_plots)
% dfile from 3dvolreg -dfile dfile
if ~exist('separate_plots','var')
    separate_plots=false;
end


if isa(dfile,'char')
    d=load(dfile);
    d=d(:,2:7);
else
    d= dfile;
end
%lbl={'Roll (IS)','Pitch (RL)','yaw (AP)','IS','RL','AP'};
% the
lbl={'IS','RL','AP','IS','RL','AP'};
% Test: the new image shifted to I, L, P
% Motion parameters were: S; R, and A; so the motion parameters are the amount of motion from new to base.
%
clr=lines(3);
figure;
if ~separate_plots
    subplot(1,2,1);
    for i=1:3
        plot(d(:,3+i),'.-','Color',clr(i,:));
        hold on;
        if i==3
            xlabel('Repetition');
        end
        
    end
    ylabel('Translation (mm)');
    xlim([0.5,size(d,1)+0.5]);
    title('Translation');
    legend(lbl(1:3));
    set(gca,'FontSize',12);
    subplot(1,2,2);
    for i=1:3
        plot(d(:,i),'.-','Color',clr(i,:));
        hold on;
        if i==3
            xlabel('Repetition');
        end
        
    end
    
    ylabel('Rotation (deg)');
    xlim([0.5,0.5+size(d,1)]);
    title('Rotation');
    legend(lbl(1:3));
    
    set(gca,'FontSize',12);
else
    for i=1:3
        subplot(3,2,i*2-1);
        plot(d(:,3+i),'.-');
        
        ylabel(['Tra (mm) - ',lbl{3+i}]);
        %  xlabel('Image #');
        
        % title(lbl(3+i));
        xlim([0.5,size(d,1)+0.5]);
        if i==1
            title('Translation');
        end
        if i==3
            xlabel('Repetition');
        end
    end
    
    for i=1:3
        subplot(3,2,i*2);
        plot(d(:,i),'.-');
        
        ylabel(['Rot (deg) - ',lbl{i}]);
        %   xlabel('Image #');
        
        %    title(lbl(i));
        xlim([0.5,0.5+size(d,1)]);
        if i==1
            title('Rotation');
        end
        if i==3
            xlabel('Repetition');
        end
    end
end
%set(gcf,'Position',[192   223   681   259]);

if isa(dfile,'char')
    prefix=strtok(dfile,'.');
    savetiff(prefix);
end
