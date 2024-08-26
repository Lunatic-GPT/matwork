function res2=merge_clusters_across_slices(m)
% merge clusters across slices, label clusters within the same PVS with the
% same number


res=m*0;

center=cell(1,size(m,3));  % stores center locations for the clusters

npvs=0;

thr=5;  % the maximum in-plane distance the clusters can have in order to belong to the same PVS 

% find the center of clusters in each slice

npvsa=zeros(1,size(m,3));

for i=1:size(m,3)
    

    tmp=clusterize2(m(:,:,i),1);

    res(:,:,i)=tmp;
    npvsa(i)=max(tmp(:));
   
    for j=1:npvsa(i)
        
        ind_tmp=ind2subb(size(tmp),find(tmp==j));
        center{i}(j,:)=mean(ind_tmp,1);
    end
    
end

if length(unique(npvsa))~=2
    error('Number of clusters different for different slices');
end

npvs=max(npvsa);

pvs_clusterID = zeros(npvs,size(m,3));  % this array will grow



% connect clusters to its nearest neighbor
for i=1:size(m,3)-1
    
    % calcuate the distance pairs
    
 
    if npvsa(i)==0
        continue;
    end
    
    
    if i==1 || sum(pvs_clusterID(:,i-1),1)==0  % start point of the PVS
        pvs_clusterID(:,i)=1:npvs;
        continue;
    end
    
    
    cind_nn=zeros(1,npvs);
    for j=1:npvs
        
        cid=pvs_clusterID(j,i-1);
        
        dist=sos(repmat(center{i-1}(cid,:),[npvs,1])-center{i},2);
        [tmp,cind_nn(j)]=min(dist);
        
    end
    
    cind1_amb=[];
    cind2_amb=[];
    
    for j=1:npvs
        
        
        if sum(cind_nn==cind_nn(j))==1  % no ambiguity
            pvs_clusterID(j,i)=cind_nn(j);
        else
            cind1_amb(end+1)=pvs_clusterID(j,i-1);
            cind2_amb(end+1)=cind_nn(j);
        end
    end
    
       %% solve the ambiguous case with brute force
        
        cind2_amb=[unique(cind2_amb),setdiff(1:npvs,cind_nn)];
        
        if length(cind2_amb)~=length(cind1_amb)
            error('ambiguous clusters not the same between two neighboring slices');
        end
        
        per=perms(1:length(cind2_amb));
        
    dist=per.*0;
    for j=1:size(per,1)
        for k=1:size(per,2)
            dist(j,k)=sos(center{i-1}(cind1_amb(k),:)-center{i}(cind2_amb(per(j,k)),:),2); 
        end
    end
    sdist=sum(dist,2);
    [tmp,ind_min]=min(sdist);
    
    for k=1:size(per,2)
        
        tmp=pvs_clusterID(:,i-1);
        pvs_clusterID(tmp==cind1_amb(k),i)=cind2_amb(per(ind_min,k));
        
    end
end

res2=res.*0;
for i=1:size(m,3)
        res2_tmp=res2(:,:,i);
    
    for j=1:npvs
        
        if pvs_clusterID(j,i)==0
            continue;
        end
        res2_tmp(res(:,:,i)==pvs_clusterID(j,i))=j;
        
        
    end
    res2(:,:,i)=res2_tmp;
end

fprintf('%d pvs found\n',npvs);








