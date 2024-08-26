function ph=phase2d(d,dim)

  ph=zeros(size(d));
  
  if dim==1
      for i=1:size(d,2)
          ph(:,i)=phase(d(:,i));
      end
  elseif dim==2
      for i=1:size(d,1)
          ph(i,:)=phase(d(i,:));
      end
  end