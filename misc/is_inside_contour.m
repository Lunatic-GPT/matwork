function res=is_inside_contour(pos,roi_pos)
        %% test whether pos is inside the contour defined by roi_pos in a plane.
        
        ph=zeros(1,size(roi_pos,1));
        
  for i=1:size(roi_pos,1)
      
      ph(i)=atan2(roi_pos(i,2)-pos(2),roi_pos(i,1)-pos(1));
      
  end
  
  dph=diff(ph);
  for i=1:length(dph)-1
      
    if sign(mod_phase(dph(i)))~=sign(mod_phase(dph(i+1)))
        res=false;
        return;
    end
        
  end
  
  res=true;
  
  
      
      function ph=mod_phase(ph)
      
      
      ph=mod(ph,2*pi);
      
      if ph>pi
          ph=ph-2*pi;
      end