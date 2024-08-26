function val=mean_roi_multVal(data,mask)
% val=mean_roi_multVal(data,mask)
% do average over roi defined in mask.
% data can be file name or 4d matrix
% mask can be file name or 3d matrix

    if isa(data,'char')
        data = ri(data);
    end
    
    if isa(mask,'char')
        mask = ri(mask);
    end

    nd=ndims(data);
    nm=ndims(mask);
    lm=length(mask(:));
    ld=length(data(:));
    
    data2=reshape(data,[lm,ld/lm]);
    
    val=zeros(max(mask(:)),ld/lm);
    for i=1:max(mask(:))
      val(i,:)=mean(data2(mask(:)==i,:),1);
    end
    
    
    sz=size(data);
    rshpsz=[max(mask(:)),sz(nm+1:end)];
    
    if length(rshpsz)>1
     val=reshape(val,rshpsz);
    end
    
    
    
    