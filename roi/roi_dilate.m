function res = roi_dilate(points,sz,rad)
% dilate mask by adding up to the rad'th neighbors
% points: n*3
res=points;
for i=-rad:rad
    for j=-rad:rad
        for k=-rad:rad
            
            pos2=points+[i,j,k];
            
            
            ind= pos2(:,1)<1 | pos2(:,1)>sz(1) |pos2(:,2)<1 | pos2(:,2)>sz(2)|pos2(:,3)<1 | pos2(:,3)>sz(3);
            
            pos2(ind,:)=[];
            
           res=cat(1,res,pos2);    
            
            
        end
    end
end

res=unique(res,'rows');
