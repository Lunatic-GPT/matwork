    function z2=permute34(z)
        sz=size(z);
        z2=zeros([sz(1:2),sz(4),sz(3)]);
        
        for i=1:sz(4)
            disp(i);
          z2(:,:,i,:)=z(:,:,:,i);
        end