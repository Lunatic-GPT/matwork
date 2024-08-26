function img2=downsample_image(img,newsize)
%downsample_image(img,newsize)

dx=1./size(img);
dx2=1./newsize;
r=dx2./dx;
if length(newsize)==2
    newsize(3)=1;
    r(3)=1;
end

img2=zeros(newsize);

for i=1:newsize(1)
  for j=1:newsize(2)
      for k=1:newsize(3)
    
          i2=round(r.*[i,j,k]);
          i1=round(r.*[i-1,j-1,k-1])+1;
          
          
          d1=img(i1(1):i2(1),i1(2):i2(2),i1(3):i2(3));
          
          img2(i,j,k)=mean(d1(:));
          
      end
  end
end

