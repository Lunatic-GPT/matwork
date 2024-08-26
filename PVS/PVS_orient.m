function pvs=PVS_orient(pvs,orient)
% dir_seg: in physical dimension (scaled by voxel size), x, y, z along the first,
% second, and third data dimensions;
voxsize=orient.voxsize;
rotmat=orient.rotmat;
npath=length(pvs);

for i=1:npath
  
    subind=pvs(i).subind;  

    posc=rotmat*subind'.*voxsize(:);
    
    orient=pca(posc');
    
    
    pvs(i).norm=orient(:,1);
    
     
    pos=posc(:,1:pvs(i).nvox_path);

    pvs(i).dir_seg=NaN(3,size(pos,2));
    
    for j=1:size(pos,2)
        
        if size(pos,2)>=5
            if j-2>0 && j+2<size(pos,2)
                orient=pca(pos(:,j-2:j+2)');
            elseif j<=2
                orient=pca(pos(:,1:5)');
            else
                orient=pca(pos(:,end-4:end)');
            end
                
            pvs(i).dir_seg(:,j)=orient(:,1);
       
        end
        
    end

 
end


