function [val,eval,val_clust]=mean_roi_clusters(data,mask)
% [val,nv]=mean_roi(data,mask)
% do average over all clusters in mask.
% data can be file name or 4d matrix
% mask can be file name or 3d matrix

tic;


    if isa(data,'char')
        data = ri(data);
    end
    
    if isa(mask,'char')
        mask = ri(mask);
    end

    m2= clusterize2(mask,1);
    
    nc=max(m2(:));
    
    val_clust=[];
    for i=1:nc
      val_clust(i,:)=mean_roi(data,m2==i);      
    end
    
    
    val=mean(val_clust,1);
    eval=std(val_clust,[],1);
    
%disp([mfilename ' finish in ', num2str(toc), ' s']);