function [pvs,mout,mout_path]=PVSPath(m,nvox_path_thr)
% m: pvs mask
% nvox_path: the minimun number of voxels in the path, clusters below the
% threshold will be removed from the output mask.
% [pvs,mout]=PVSPath(m,nvox_path)
% pvs: a struct array 
% mout: the mask after thresholding
% mout_path: the mask containing only the path voxels

[~,subind_c]=clusterize2(m>0,2);
tic;
mout=m*0;
mout_path=m*0;
n_rm=0;
n_kp=0;
for i=1:length(subind_c)
    [nvox_path,subind]=calcPath(subind_c{i});
    
    if nvox_path>=nvox_path_thr
      n_kp=n_kp+1;
      pvs(n_kp).nvox_path=nvox_path;
      pvs(n_kp).subind=subind;
      pvs(n_kp).ind=sub2indb(size(m),subind);
      mout(pvs(n_kp).ind)=1;
      mout_path(pvs(n_kp).ind(1:nvox_path))=1;
      
    else
      n_rm=n_rm+1;
    end
end


fprintf('%d (%d) PVS found (removed)\n',length(subind_c),n_rm);


function [nvox_path,pos,maxlen]=calcPath(pos)


a=pca(pos);
v=a(:,1);

d=zeros(1,size(pos,1));

for i=1:size(pos,1)
  d(i)=pos(i,:)*v;
end

[tmp,imax]=max(d);
[tmp,imin]=min(d);

pair=[imin,imax];

[s,w]=find_neighbors2(pos);

sp=sparse(s(:,1)',s(:,2)',w);
sp(size(sp,1)+1:size(sp,2),:)=0;
sp=tril(sp + sp');
[maxlen,pathmax] = graphshortestpath(sp,pair(1),pair(2),'Directed',false);
nvox_path=length(pathmax);
xpath=setdiff(1:size(pos,1),pathmax);
pos=[pos(pathmax,:);pos(xpath,:)];



function [p,w] = find_neighbors2(pos)

p=[];
w=[];
for i=1:size(pos,1)
    pos1=pos(i,:);
    for j=i+1:size(pos,1)
        pos2=pos(j,:);
        if ~any(abs(pos1-pos2)>1)
            p(end+1,:)=[i,j];
            w(end+1)= sqrt(sum(abs(pos1-pos2).^2));
        end
    end
end


  

function subind=ind2subb(mtrx,ind_arr)

subind=zeros(length(ind_arr),length(mtrx));
for j=1:length(ind_arr)
    ind=ind_arr(j);
    for i=1:length(mtrx)
        
        if i==length(mtrx)
            subind(j,1)=ind;
        else
            subind(j,end-i+1)=floor((ind-1)/prod(mtrx(1:end-i)))+1;
        end
        ind=ind-(subind(j,end-i+1)-1)*prod(mtrx(1:end-i));
    end
    
end







function [omask,clusters]=clusterize2(mask,thr)
%[mask,clusters]=clusterize2(mask,thr)
% find connected clusters in the mask and remove voxels with cluster size
% less than thr.
% only version; NO clusterize or clusterize1.


if ~exist('thr','var')
    thr = 1;
end

omask = zeros(size(mask));
sz = size(mask);

if length(sz)==2
    sz(3)=1;
end

clusters = {};
cluster_temp = zeros(sz(1)*sz(2)*sz(3),3);


ind=find(mask);
iclust=1;
for i=1:length(ind)
    ii=ind2subb(size(mask),ind(i));
    if length(ii)==2
        ii(3)=1;
    end
    if(mask(ii(1),ii(2),ii(3))==0)
        continue;  % already in a cluster.
    end
    
    cluster_temp(1,:)=ii;
    
    
    mask(ii(1),ii(2),ii(3))=0;
    nvc=1; % number of voxels in a cluster.
    
    ivc = 1;
    while ivc<=nvc
        
        nb = find_neighbors_clust(cluster_temp(ivc,:),mask);
        
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

nvox=zeros(1,length(clusters));

for i=1:length(clusters)
    nvox(i)=size(clusters{i},1);
end

[nvox_sort,ind]=sort(nvox,'descend');

clusters=clusters(ind);

clusters=clusters(nvox_sort>=thr);

nclus = length(clusters);

for i=1:nclus
    
    for j=1:size(clusters{i},1)
        
        a=clusters{i}(j,:);
        omask(a(1),a(2),a(3)) =i;
        
        
        
    end
    
end
% sort clusters to be compatible with earlier version

    for j=1:length(clusters)
        ind=sub2indb(size(omask),clusters{j});
        
        [tmp,i_ind]=sort(ind);
        clusters{j}=clusters{j}(i_ind,:);
    end



%  fprintf('number of voxels in mask = %d\n',length(find(omask)));



function nb = find_neighbors_clust(pos,mask)

nb=[];
sz = size(mask);
if numel(sz) <3
    sz(3)=1;
end

for i=-1:1
    for j=-1:1
        for k=-1:1
            if(abs(i)==1 &&abs(j)==1 && abs(k)==1)
                %  continue;
            end
            
            pos2=pos+[i,j,k];
            if pos2(1)<1 || pos2(1)>sz(1) ||pos2(2)<1 || pos2(2)>sz(2)||pos2(3)<1 || pos2(3)>sz(3)
                continue;
            end
            
            if mask(pos2(1),pos2(2),pos2(3))>0
                nb(end+1,:)=pos2;
            end
        end
    end
end













