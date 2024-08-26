function res=find_mask_near_path(path,data,thr,rad)
%res=find_mask_near_path(path,data,thr,rad)
% path: n*3
% Find voxels data>thr or data<abs(thr) (if thr<0) within rad voxels
% surrounding the path
% If more than one cluster is found, only the one connected with path will
% be returned;
% The voxels in path will always be included, even if they are below the
% threshold.
% output: m*3
path_dil=roi_dilate(path,size(data),rad);
ind_new=sub2indb(size(data),path_dil);
if thr>0
res=path_dil(data(ind_new)>thr,:);
else
    res=path_dil(data(ind_new)<abs(thr),:);
end
res=unique(cat(1,res,path),'rows');

res=keep_connected(res,path(1,:));


function res=keep_connected(mask,seed)
   
%only keep the cluster connected with seed
  
a=min(mask,[],1);
b=max(mask,[],1);
sz=b-a+1;
m=zeros(sz);
m=setv_points(m,mask-a+1,1);

[~,~,clust]=clusterize2(m);
for i=1:length(clust)
   if any(sum(clust{i}==seed-a+1,2)==3)
       
       res=clust{i}+a-1;
       break;
   end
       
    
end
