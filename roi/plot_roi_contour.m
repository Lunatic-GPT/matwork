function plot_roi_contour(roi,color,linewidth)
% plot_roi_contour(roi,color,linewidth)

hold on;
for i=2:size(roi,1)-1
    for j=2:size(roi,2)-1
        if roi(i,j)>0  && roi(i,j-1)==0  % left edge of a pixel
            plot(j+[-0.5,-0.5],i+[-0.5,0.5],'Color',color,'LineWidth',linewidth);
        end
        
        if roi(i,j)>0  && roi(i,j+1)==0  % right edge of a pixel
            plot(j+[0.5,0.5],i+[-0.5,0.5],'Color',color,'LineWidth',linewidth);
        end
        
        if roi(i,j)>0  && roi(i-1,j)==0  % top edge of a pixel
            plot(j+[-0.5,0.5],i+[-0.5,-0.5],'Color',color,'LineWidth',linewidth);
        end
        
        if roi(i,j)>0  && roi(i+1,j)==0  % bottom edge of a pixel
            plot(j+[-0.5,0.5],i+[0.5,0.5],'Color',color,'LineWidth',linewidth);
        end
    end
end
            

