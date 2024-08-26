function [rows,cols] = getCalibSizeTSE(mask,szmin)
% 

imin=0;
imax=size(mask,1);
jmin=0;
jmax=size(mask,2);

found=false;
for i=1:size(mask,1)-szmin(1)+1
    for j=1:size(mask,2)-szmin(2)+1
        
      if sum(vec(mask(i:i+szmin(1)-1,j:j+szmin(2)-1)))==prod(szmin)
          if ~found
              imin=i;
              jmin=j;
              found=true;
          else
              imax=i;
              jmax=j;
          end
              
      end
      
    end
end

rows=imin:imax+szmin(1)-1;
cols=jmin:jmax+szmin(2)-1;



