function [res,res_clust]=sen_ppv_dsc(m,gt)
%  [res,res_clust]=sen_ppv_dsc(m,gt)
% res and res_clust both contain sen, ppv, and dsc
%res_clust for cluster level sen, ppv, and dsc
        
    tp=sum(vec(gt>0&m>0));
    fp=sum(vec(gt==0&m>0));
    fn=sum(vec(gt>0&m==0));
    
    res.sen=tp/(tp+fn);
    res.ppv=tp/(tp+fp);
    res.dsc=2*tp/(2*tp+fp+fn);
    res.nvox=sum(m(:)>0);
    res.nvox_gt=sum(gt(:)>0);
     
    if nargout==1
        return;
    end
    %%
    n_tp=0; % number of true positive clusters
    n_fp=0; % number of false positive clusters
    n_fn=0;
    
    c=clusterize2(m);
    c_gt=clusterize2(gt);
    
    m_fp=c*0;
    tic;
    for j=1:max(c(:))
        m=(c==j)&gt>0;
        if any(m(:))
            n_tp=n_tp+1;
        else
            m_fp(c==j)=1;
           n_fp=n_fp+1; 
        end
        time_left(j,max(c(:)),toc);
    end
      tic;
    m_fn=c_gt*0;
    for j=1:max(c_gt(:))
        m=(c_gt==j)&c>0;  
        if ~any(m(:))
           n_fn=n_fn+1;
           m_fn(c_gt==j)=1;
        end   
        time_left(j,max(c(:)),toc);
    end
      
    res_clust.sen=n_tp/(n_tp+n_fn);
    res_clust.ppv=n_tp/(n_tp+n_fp);
    res_clust.dsc=2*n_tp/(2*n_tp+n_fp+n_fn);
   res_clust.nclust_gt=max(c_gt(:));
   res_clust.nclust=max(c(:));
   
    
