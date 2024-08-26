function process_pvs(fname,skip_q)
if ~exist('skip_q','var')
    skip_q=false;
end

fname=name4pat(fname);
for i=1:length(fname)
    process_pvs_1file(fname{i},skip_q);    
end

function process_pvs_1file(fname,skip_q)
% 6/10/2019: in PVSPath, delete ind (n*1); use (indc n*3) instead to
% speed up. delete path variable, save nvox_path; since
% path{i}=1:nvox_path(i)


if ~exist('rmBranch','var')
    rmBranch=false;
end


[dname,fname2]=fileparts(fname);
prefix=remove_suffix(fname2,'.nii');

%prefix=strtok(fname2,'.');
if ~exist(prefix,'dir')
    mkdir(prefix);
end

nii=load_untouch_niigz(fname);

m=nii.img>0.8;
sz=size(m);
cd(prefix);
prefix2=[filename(prefix),'_PVS'];
pvs=PVSPath(m);
save_mask(pvs,sz,nii,prefix2);
 

if skip_q
    cd ..;
    return;
end
  
orient=get_orient_from_nii(nii);

pvs=PVS_orient(pvs,orient);
PVS_quantify(pvs,nii,filename(prefix));



cd('..');
%%

function save_mask(pvs,sz,nii,prefix)

[c,c_path]=set_c_from_ind(pvs,sz);
nii.img=c;
save_untouch_niigz(nii,prefix);

prefix2=[prefix,'Path'];
nii.img=c_path;
save_untouch_niigz(nii,prefix2);


function pvs=PVSPath(m)
% m: pvs mask
% prefix: file name of the output file
global voxsize;

[~,subind_c]=clusterize2(m>0,2);
tic;
for i=1:length(subind_c)
    [nvox_path,subind]=calcPath(subind_c{i});
    pvs(i).nvox_path=nvox_path;
    pvs(i).subind=subind;
    pvs(i).ind=sub2indb(size(m),subind);
end
fprintf('%d PVS found\n',length(pvs));


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



function [c,c_path]=set_c_from_ind(path,sz)
        
c=zeros(sz,'int16');
c_path=zeros(sz,'int16');

for i=1:length(path)    
    c(path(i).ind)=i;    
    c_path(path(i).ind(1:path(i).nvox_path))=i;
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


function pos2=connect_pos(pos)

pos2=pos(1,:);
    
    n=pos(2,:)-pos(1,:);  
    npix=ceil(sqrt(sum(n.^2)))*2;
    
    for j=1:npix-1
        d=n/npix*j;
        
        if ~any(abs(d)-0.5<0.001)
            pos2(end+1,:) = round(pos(1,:)+d);
        end
    end
    
    pos2=[pos2;pos(2,:)];
    pos2=unique(pos2,'rows','stable');







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













