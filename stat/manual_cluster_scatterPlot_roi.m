function m_out=manual_cluster_scatterPlot_roi(x,y,mask)
% mask contains values from 0 to nroi


h=impoly(gca,[]);

ButtonName = questdlg('ROI OK', 'Confirm ROI','Yes', 'No','Yes');

if strcmp(ButtonName,'No')
    return;
end

api = iptgetapi(h);
roi_pos = api.getPosition();
m_out=0*mask;
nroi=max(mask(:));

for i=1:nroi
    
    if is_inside_contour([x(i),y(i)],roi_pos)
        m_out(mask==i)=i;
    end
end

save manual_cluster_scatterPlot_roi m_out

  
  
  
      
      
      
      
      
  
  
      
  
  
  