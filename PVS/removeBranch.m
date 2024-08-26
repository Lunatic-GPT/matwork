function [pvs,c_hasBranch,n_hasBranch]=removeBranch(pvs,sz)


%% if there is a gap between the voxel and its path, delete that voxel

c_hasBranch=zeros(sz,'int16'); %1: voxels that were removed; 2: other voxels from the same cluster
c_removed=zeros(sz,'int16');
n_hasBranch=0;
npvs=length(pvs);
for i=1:length(pvs)
  voxremoved=false;
    while 1
       
        voxremoved_thisitr=false;
        nvox=pvs(i).nvox_path;
        while  nvox<length(pvs(i).ind)
            
            dist=sos(pvs(i).subind(nvox+1,:)-pvs(i).subind(1:pvs(i).nvox_path,:),2);
            if min(dist)>=2
                
                [tmp,indtmp]=min(dist);
               
                if hasgap(pvs(i).subind,nvox+1,indtmp)                    
                    if ~voxremoved  %only assign to 2 and increment npvs_hasBranch once
                        val=c_hasBranch(pvs(i).ind);
                        c_hasBranch(pvs(i).ind(val~=1))=2;  %1: branch voxel; 2: other voxels in the PVS
                        voxremoved=true;
                        n_hasBranch=n_hasBranch+1;
                    end
                    
                    c_hasBranch(pvs(i).ind(nvox+1))=1;
                    
                    c_removed(pvs(i).ind(nvox+1))=1;
                %    c(pvs(i).ind(nvox+1))=0;
                    
                    pvs(i).ind(nvox+1)=[];
                    pvs(i).subind(nvox+1,:)=[];        
                    voxremoved_thisitr=true;
                    continue;
                end
            end
            nvox=nvox+1;
        end
        
        if ~voxremoved_thisitr
            break;
        end
        %if voxremoved, check again.
    end
   
end

fprintf('%d out of %d PVSs have branches\n',  n_hasBranch,npvs);
    
% add the removed voxels back as new PVSs
%{
if n_hasBranch>0
    pvs2=PVSPath(c_removed);
    [pvs2,c_hasBranch2,n_hasBranch2]=removeBranch(pvs2,sz);
    
    c_hasBranch(c_hasBranch2>0)=c_hasBranch2(c_hasBranch2>0);
    pvs=[pvs,pvs2];
    n_hasBranch=n_hasBranch+n_hasBranch2;
end
%}

function res=hasgap(pos,i,j)

pos2=connect_pos(pos([i,j],:));

for i=2:size(pos2,1)-1
    
    if ~any(pos(:,1)==pos2(i,1) &pos(:,2)==pos2(i,2) & pos(:,3)==pos2(i,3))
        res= true;
        return;
    end
end
res=false;



function pos2=connect_pos(pos)

pos2=pos(1,:);

 n=pos(2,:)-pos(1,:);
 
 npix=ceil(sqrt(sum(n.^2)))*2;
 
 for j=1:npix-1
    
     dpos=j*n/npix;
     
     if ~any(abs(abs(dpos)-0.5)<0.00001)
     
         tmp=round(pos(1,:)+j*n/npix);
         pos2(end+1,:) = tmp; 
         
     end
    
    
 end
 
 pos2(end+1,:)=pos(2,:);
 
 pos2=unique(pos2,'rows','stable');
 
    
