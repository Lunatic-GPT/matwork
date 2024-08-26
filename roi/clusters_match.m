function [iab,iab_ol,iab_mindist]=clusters_match(a,b,max_clust_dist,do_clusterize)
% [iab,iab_ol]=clusters_match(a,b,max_clust_dist,do_clusterize)
% a and b: masks
% do_clusterize: whether to re-label the clusters; needed if the non-zero
% voxels all have the same value
% max_clust_dist: maximum distance between centers of overlapping clusters;
% (unit: pixels)
%iab: (n*2); pairs of index values that match in space;
%iab_ol: pairs of index values that match in space but one of them match with multiple ROIs on the other mask;
if ~exist('do_clusterize','var')
    do_clusterize=false;
end

%if max(a(:))==1 || max(b(:))==1
if do_clusterize
    disp('clusterize the rois in clusters_match');
    a=clusterize2(a);
    b=clusterize2(b);
end

a(isnan(a))=0;
b(isnan(b))=0;

%end

for i=1:max(a(:))
    if sum(a(:)==i)>0
        com_a(i,:)=mean(roiCOM(a==i),1);
    end
end

for i=1:max(b(:))
    if sum(b(:)==i)>0
        com_b(i,:)=mean(roiCOM(b==i),1);
    end
end

b2=b*0;
a2=a*0;
nmatch=0;
ia=[];
ib=[];
dist=[];
for i=1:max(a(:))
    
    nfound=0;
    for j=1:max(b(:))
        if ~any(com_a(i,:)>0) || ~any(com_b(j,:)>0)
            continue;
        end
        
        if sos(com_a(i,:) - com_b(j,:)) <=max_clust_dist
            ia(end+1)=i;
            ib(end+1)=j;
            
            dist(end+1)=sos(com_a(i,:) - com_b(j,:));
            nfound=nfound+1;
        end
    end
    
end
ia=ia(:);
ib=ib(:);
[iab,ind_rm]=remove_ambiguity(ia,ib);
iab_ol=[ia(ind_rm),ib(ind_rm)];

iind_md=resolve_ambiguity(ia(ind_rm),ib(ind_rm),dist(ind_rm));

iab_mindist=[ia(ind_rm(iind_md)),ib(ind_rm(iind_md))];

%[iab,
fprintf('Total clusters = %d (a) %d (b);\n',...
    length(unique(a(:)))-1,length(unique(b(:)))-1);
fprintf('Match non-ambiguous: %d \n',size(iab,1));
fprintf('Match ambiguous: %d (mask a); %d (mask b); \n',...
    length(unique(iab_ol(:,1))),length(unique(iab_ol(:,2))));
fprintf('Resolved ambiguous: %d \n',...
    size(iab_mindist,1));

function  [ind_mindist,ind_rm]=resolve_ambiguity(ia,ib,dist)

ind_mindist=[];
ind_rm=[];
while length(ind_rm)+length(ind_mindist)<length(ia)
    
    
    ind_clust=ind_first_connected_pairs(ia,ib,[ind_mindist;ind_rm]);
    
    [~,ii]=min(dist(ind_clust));
    ind_mindist=[ind_mindist;ind_clust(ii)];
    
    ind2rm=find(ia==ia(ind_clust(ii)) | ib==ib(ind_clust(ii)));
    ind2rm(ind2rm==ind_clust(ii))=[];
    ind_rm=[ind_rm;ind2rm];
    
end



function ind=ind_first_connected_pairs(ia,ib,ind_exc)

ind_new=1:length(ia);
ind_new(ind_exc)=[];
ia(ind_exc)=[];
ib(ind_exc)=[];

iind=1;
  while 1
           
      ind_tmp1=find_ind_same(ia(iind),ia);
      ind_tmp2=find_ind_same(ib(iind),ib);
      
      if length(ind_tmp1)==length(iind) && length(ind_tmp2)==length(iind)
          break;
      else
          iind=unique([ind_tmp1;ind_tmp2]);
      end
  end

 ind=ind_new(iind);
  
function ind=find_ind_same(sub,group)
        %find all indices for values in group that are contained in sub
        ind=[];
        for i=1:length(sub)
            tmp=find(group==sub(i));
           ind=[ind;tmp]; 
        end
        ind=unique(ind);
        
        
        

function [iab,ind_rm]=remove_ambiguity(ia,ib)


if length(ia)>length(unique(ia))
    fprintf('multiple matches for some clusters in mask 1.\n');
end

if length(ib)>length(unique(ib))
    fprintf('multiple matches for some clusters in mask 2.\n');
end


ind_rm=[];
for i=1:length(ia)
    itmp=find(ia==ia(i));
    if length(itmp)>1
        % fprintf('Cluster %d removed from mask 1\n',ia(i));
        ind_rm=[ind_rm;itmp];
    end
    
    itmp=find(ib==ib(i));
    if length(itmp)>1
        % fprintf('Cluster %d removed from mask 2\n',ib(i));
        ind_rm=[ind_rm;itmp];
    end
end
ind_rm=unique(ind_rm);
if ~isempty(ind_rm)
fprintf('Ambiguous matches from mask 1:\n');
disp(ia(ind_rm));

fprintf('from mask 2:\n');
disp(ib(ind_rm));
end

iab=[ia(:),ib(:)];
iab(ind_rm,:)=[];
