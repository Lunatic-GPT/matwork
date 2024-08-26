function plot_motion_afni_mat(dfile)
% dfile is a mat file; generated after 

d=load(dfile);
%lbl={'Roll (IS)','Pitch (RL)','yaw (AP)','IS','RL','AP'};
% the 
% lbl={'IS','RL','AP','IS','RL','AP'};
% Test: the new image shifted to I, L, P
% Motion parameters were: S; R, and A; so the motion parameters are the amount of motion from new to base.
% 
% figure;
%     subplot(1,2,1);    
%     plot(d.an');
%     hold on;
%     plot(squeeze(abs(d.total(3,:))),'k-');
%     ylabel('Rotation (Deg)');
%     xlabel('Image #');
%     legend(lbl(1:3));
%     
%     xlim([1,size(d.an,2)]);
%     subplot(1,2,2);    
%     plot(squeeze(d.shift([3,1,2],1,:))');
%     
%     ylabel('Translation (mm)');
%     xlabel('Image #');
%     legend(lbl(4:6));
%     xlim([1,size(d.an,2)]);
%     
%     set(gcf,'Position',[192   223   681   259]);
    
    prefix=strtok(dfile,'.');
    d=cat(2,d.an',squeeze(d.shift([3,1,2],1,:))');
    do_plot(d,prefix,false);
    
function do_plot(d,prefix,separate_plots)
        
lbl={'IS','RL','AP','IS','RL','AP'};
% Test: the new image shifted to I, L, P
% Motion parameters were: S; R, and A; so the motion parameters are the amount of motion from new to base.
%
clr=lines(3);
figure;
if ~separate_plots
     subplot(2,1,1);
     for i=1:3    
      plot(d(:,3+i),'.-','Color',clr(i,:));
      hold on;
     end
      set(gca,'FontSize',12);
     ylabel('(mm)');
     xlim([0.5,size(d,1)+0.5]);
     title('Translation');
     legend(lbl(1:3));
    set(gca,'xtick',0:20:160);
     set(gca,'ytick',-1:0.2:1);
     subplot(2,1,2);
    for i=1:3
      plot(d(:,i),'.-','Color',clr(i,:));
      hold on;
    end
     set(gca,'FontSize',12);
     ylabel('(degree)');
     xlim([0.5,0.5+size(d,1)]);
     title('Rotation');
     legend(lbl(1:3));
     xlabel('TR index');
    set(gca,'xtick',0:20:160);
     set(gca,'ytick',-1:0.2:1);
    
    
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
        
    end
end
%set(gcf,'Position',[192   223   681   259]);
set(gcf,'Position',[10.5521    3.2604    3.5625    4.9063]);
savetiffc(prefix);
