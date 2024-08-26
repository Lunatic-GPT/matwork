function a=shift_image_fracvox(a,shift)
%a=shift_image_fracvox(a,shift)
% shift: 1*3 vector.  positive to shift to larger indices and negative to
% smaller indices.  Same convention as in circshift.
    

if shift(1)~=0 
    if  floor(shift(1))==ceil(shift(1))
    
      a=circshift(a,[ceil(shift(1)),0,0,0]);
    else
        
        
      a1=circshift(a,[ceil(shift(1)),0,0,0]);
      a2=circshift(a,[floor(shift(1)),0,0,0]);
    
      w2=ceil(shift(1))-shift(1);
      w1=shift(1)-floor(shift(1));
      a=a1*w1+a2*w2;
      
    end
end


%if shift(2)~=0 
   %{ 
if  floor(shift(2))==ceil(shift(2))
    
      a=circshift(a,[0,ceil(shift(2)),0,0]);
    else
     %}   
      %{  
      a1=circshift(a,[0,ceil(shift(2)),0,0]);
      a2=circshift(a,[0,floor(shift(2)),0,0]);
    
      w2=ceil(shift(2))-shift(2);
      w1=shift(2)-floor(shift(2));
      a=a1*w1+a2*w2;
        %}
        ph=exp(1i*shift(2)*2*pi/size(a,2)*(-size(a,2)/2:size(a,2)/2-1));
        fa=fft1c(a,2);
        fa=fa.*repmat(ph,[size(a,1),1,size(a,3)]);
        a=ifft1c(fa,2);
        
        
%    end
%end


if shift(3)~=0 
    if  floor(shift(3))==ceil(shift(3))
    
      a=circshift(a,[0,0,ceil(shift(3)),0]);
    else
      a1=circshift(a,[0,0,ceil(shift(3)),0]);
      a2=circshift(a,[0,0,floor(shift(3)),0]);
    
      w2=ceil(shift(3))-shift(3);
      w1=shift(3)-floor(shift(3));
      a=a1*w1+a2*w2;
    end
end
