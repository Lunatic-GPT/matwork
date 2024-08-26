function res = mtimes(a,b)

if isa(a,'Wavelet_ISU') == 0
    error('In  A.*B only A can be Wavelet operator');
end

res=zeros(size(b));
if a.adjoint
    
    for i=1:size(b,3)
        for j=1:size(b,4)
         res(:,:,i,j) = a.wavl'*b(:,:,i,j)*a.wavr;
        end
    end
    
    if numel(a.N)==3
        for j=1:size(b,3)
        tmp=0;
          for i=1:size(b,3)
           tmp=tmp+res(:,:,i,:)*a.wavs(i,j);
          end
          res2(:,:,j,:)=tmp;
        end
        res=res2;
    end
else
    for i=1:size(b,3)
        for j=1:size(b,4)
         res(:,:,i,j) = a.wavl*b(:,:,i,j)*a.wavr';
        end
    end
    
    if numel(a.N)==3
        for j=1:size(b,3)
        tmp=0;
          for i=1:size(b,3)
           tmp=tmp+res(:,:,i,:)*a.wavs(j,i);
          end
          res2(:,:,j,:)=tmp;
        end
        res=res2;
    end
    
end


