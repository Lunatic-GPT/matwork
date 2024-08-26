function res=sen_ppv_dsc_path(m,gt,max_dist)
%  res=sen_ppv_dsc(m,gt)
% res contain sen, ppv, and dsc
% max_dist: the maximum dist between m and gt voxels.  If a neighboring
% voxel is not found within max_dist, then it is considered no match.

if ~exist('max_dist','var')
    max_dist=1;
end


    %%
    n_tp=0;
    n_fp=0;
    n_fn=0;
    
    ind_m=find(m>0);
    ind_gt=find(gt>0);
    ind_match_m=[];
    ind_match_gt=[];
    
    for i=1:length(ind_m(:))
        ind_tmp=find_neighbor(ind_m(i),gt,max_dist);
        
        if ind_tmp>0 %&&~any(ind_tmp==ind_match_m)
            n_tp=n_tp+1;
            ind_match_m(end+1,:)=[ind_m(i),ind_tmp];
        else
            n_fp=n_fp+1; 
        end
            
    end
    
   for i=1:length(ind_gt(:))
        ind_tmp=find_neighbor(ind_gt(i),m,max_dist);
        
        if ind_tmp>0 %&&~any(ind_tmp==ind_match_gt)
            ind_match_gt(end+1,:)=[ind_gt(i),ind_tmp];
        else
            n_fn=n_fn+1; 
        end
    end
    
      
    res.sen=n_tp/(n_tp+n_fn);
    res.ppv=n_tp/(n_tp+n_fp);
    res.dsc=2*n_tp/(2*n_tp+n_fp+n_fn);
    res.ind_match_m=ind_match_m;
    res.ind_match_gt=ind_match_gt;
    
function res=find_neighbor(ind,m,max_dist)
sz=size(m);
ind=ind2subb(sz,ind);
ind2=index_within_limit([ind(:)-max_dist,ind(:)+max_dist],sz);

 dist=Inf;
 res=0;
 for i=ind2(1,1):ind2(1,2)
     for j=ind2(2,1):ind2(2,2)
         for k=ind2(3,1):ind2(3,2)
            if m(i,j,k)>0
                
                dist_tmp=sos([i,j,k]-ind);
                
                if dist_tmp<dist
                    dist=dist_tmp;
                    res=sub2indb(sz,[i,j,k]); 
                end
            end
             
         end
     end
 end








        
