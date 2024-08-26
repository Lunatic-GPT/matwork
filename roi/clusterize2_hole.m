function omask=clusterize2_hole(mask,thr)

% remove holes in mask with number of voxels less than thr; 3D holes


if ~exist('thr','var')
    thr = 1;
end

mask=(mask==0);

omask = mask;
sz = size(mask);

if length(sz)==2
    sz(3)=1;
end
    
clusters = {};
cluster_temp = zeros(sz(1)*sz(2)*sz(3),3);

for i=1:sz(1)
    for j=1:sz(2)
        for k=1:sz(3)
         if mask(i,j,k) >0      
           cluster_temp(1,:) = [i,j,k];
         else
             continue;
         end

         mask(i,j,k)=0;
         nvc=1; % number of voxels in a cluster.
            
         ivc = 1;
         while ivc<=nvc 
        
          nb = find_neighbors(cluster_temp(ivc,:),mask);
        
          nnb = size(nb,1);   
         if nnb > 0 
          cluster_temp(nvc+1:nvc+nnb,:) = nb;
          nvc=nvc+nnb;
          for a=1:nnb
             mask(nb(a,1),nb(a,2),nb(a,3))=0;
          end
         end
        ivc = ivc+1;
         end
       clusters{end+1} = cluster_temp(1:nvc,:); 

        end
    end
end

nclus = length(clusters);
clus_ind = [];
for i=1:nclus
    if size(clusters{i},1) < thr
        clus_ind(end+1) = i;
        for j=1:size(clusters{i},1)
            a=clusters{i}(j,:);
            omask(a(1),a(2),a(3)) =0;
        end
        
     end
end
    
omask=omask==0;       
    
function nb = find_neighbors(pos,mask)

   sz = size(mask);
   if numel(sz) <3
       sz(3)=1;
   end
   nnb = 0;
   nb = [];
  % nb = zeros(6,3);
   if pos(3) < sz(3) && mask(pos(1),pos(2),pos(3)+1)>0    
       nb(nnb+1,1:3)=pos+[0,0,1];
       nnb = nnb+1;
   end
   
   if pos(3) > 1 && mask(pos(1),pos(2),pos(3)-1)>0    
       nb(nnb+1,1:3)=pos-[0,0,1];
       nnb = nnb+1;
   end
   
   if pos(2) < sz(2) && mask(pos(1),pos(2)+1,pos(3))>0    
       nb(nnb+1,1:3)=pos+[0,1,0];
       nnb = nnb+1;
   end
   
   if pos(2) >1 && mask(pos(1),pos(2)-1,pos(3))>0    
       nb(nnb+1,1:3)=pos-[0,1,0];
       nnb = nnb+1;
   end
   
   if pos(1) < sz(1) && mask(pos(1)+1,pos(2),pos(3))>0    
       nb(nnb+1,1:3)=pos+[1,0,0];
       nnb = nnb+1;
   end
   
   if pos(1) >1 && mask(pos(1)-1,pos(2),pos(3))>0    
       nb(nnb+1,1:3)=pos-[1,0,0];
   end
%nb = nb(1:nnb,:);   
   
     
   
