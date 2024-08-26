function voxpos2mask(vox,pos,prefix_parent,prefix_out)
%voxpos2mask(vox,pos,prefix_parent,prefix_out)
% all voxels intersect with the roi will be selected.
[err,info]=BrikInfo([prefix_parent,'+orig']);

ad=abs(info.DELTA);
x=linspace(-vox(1)*0.49,vox(1)*0.49,ceil(vox(1)/ad(1))+1)+pos(1);
y=linspace(-vox(2)*0.49,vox(2)*0.49,ceil(vox(2)/ad(2))+1)+pos(2);
z=linspace(-vox(3)*0.49,vox(3)*0.49,ceil(vox(3)/ad(3))+1)+pos(3);

d=info.DATASET_DIMENSIONS;
mask=zeros(d(1:3));
for i=1:length(x)
    for j=1:length(y)
        for k=1:length(z)
            
            ijk=xyz2ijk([x(i),y(j),z(k)],info)+[1,1,1];
            mask(ijk(1),ijk(2),ijk(3))=1;
        end
    end
end

history = sprintf('voxpos2mask(vox=%s, pos=%s,prefix_parent=%s)',num2str(vox),num2str(pos),prefix_parent);
WriteBrikEZ(mask,info,history,prefix_out);



function ijk = xyz2ijk(xyz,info)
% ijk indices in afni.  
% This is one less than the indices in matlab.    
                dxyz = info.DELTA;
                orgn = info.ORIGIN; 
                ijk=round((xyz-orgn)./dxyz);
                
                