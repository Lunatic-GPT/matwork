function mask=clusterize2_2d_hole(mask,thr)

%mask=clusterize2_2d_hole(mask,thr)
% remove holes in mask with number of voxels less than thr; slice by slice

mask=(mask==0);
for i=1:size(mask,3)
    
    tmp=mask(:,:,i);
    if sum(tmp(:))<thr
        mask(:,:,i)=0;
        continue;
    end
    
    mask(:,:,i)=clusterize2(mask(:,:,i),thr);
end

mask=(mask==0);
