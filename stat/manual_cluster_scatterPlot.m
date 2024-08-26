function x=manual_cluster_scatterPlot(xdata,ydata,mask)


figure;plot(xdata(mask),ydata(mask),'o');
h=impoly(gca,[]);

ButtonName = questdlg('ROI OK', 'Confirm ROI','Yes', 'No','Yes');

if strcmp(ButtonName,'No')
    return;
end

api = iptgetapi(h);
roi_pos = api.getPosition();
x=0*mask;

for i=1:length(xdata(:))
    if mask(i)==0
        continue;
    end
    
    if is_inside_contour([xdata(i),ydata(i)],roi_pos)
        x(i)=1;
    else
        x(i)=2;
    end
end

save manual_cluster_scatterPlot x

  
  
  
      
      
      
      
      
  
  
      
  
  
  