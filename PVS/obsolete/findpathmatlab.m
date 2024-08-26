
%{

function [npix_path,pos2,maxlen]=findPathMatlab(pos)

%% pos is the ind of voxels in the cluster.1*nvox
%% sz: matrix size.
%%{


%pos=ind2subb(sz,ind);


for i=1:3
    [tmp,imax(i)]=max(pos(:,i));
    [tmp,imin(i)]=min(pos(:,i));
end

%iall=unique([imax,imin]);
iall=([imax,imin]);
dist=0;
% pair=[];
% for i=1:length(iall)
%     for j=i+1:length(iall)
%         
%         disttmp=sos(pos(iall(i),:)-pos(iall(j),:));
%         if dist<disttmp
%             pair=[iall(i),iall(j)];
%             dist=disttmp;
%         end
%         
%     end
% end
for i=1:6
    for j=2:6
        
        if j<=i
            continue;
        end
        
        if iall(i)==iall(j)
            continue;
        end
        
        disttmp=sos(pos(iall(i),:)-pos(iall(j),:));
        if dist<disttmp
            pair=[iall(i),iall(j)];
            dist=disttmp;
        end
        
    end
end

[s,w]=find_neighbors2(pos);
sp=sparse(s(:,1)',s(:,2)',w);


sp(size(sp,1)+1:size(sp,2),:)=0;

sp=tril(sp + sp');
[maxlen,pathmax] = graphshortestpath(sp,pair(1),pair(2),'Directed',false);


xpath=setdiff(1:size(pos,1),pathmax);

pos2=[pos(pathmax,:);pos(xpath,:)];
npix_path=length(pathmax);

%}