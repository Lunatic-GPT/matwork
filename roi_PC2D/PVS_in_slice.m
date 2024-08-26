function PVS_in_slice(f_pvs,nii_pc,prefix)
% f_pvs: the matlab file containing pvs array and orient
% f_pc: the nifti file for PC 
% prefix: the prefix of the output nii.gz file.

orient_slice=get_orient_from_nii(nii_pc);
orient_pvs=load(f_pvs);

m=PVS_in_slice_kernel(orient_pvs.pvs,orient_pvs,orient_slice);

nii_pc.img=m;
save_untouch_niigz(nii_pc,prefix);




function m=PVS_in_slice_kernel(pvs,orient_pvs,orient_slice)
% input:
% pvs: array of PVS structure
% orient_pvs: the orientation of PVS
% orient_slice: the orientation of slice
% output:
% the mask with voxels intersecting the PVS set to the indices of the PVS
% structure array.



count=0;
m=zeros(orient_slice.sz);

for i=1:length(pvs)
    sub=pvs(i).subind(1:pvs(i).nvox_path,:);
    
    xyz=ijk2xyz(sub,orient_pvs);
    
    norm=orient_slice.rotmat(:,3);
    res=dist2plane(xyz,orient_slice.center,norm);
    
    iint=find(sign(res(1:end-1)).*sign(res(2:end))<=0);
    if length(iint)>1 || isempty(iint)
        continue;
    end
    
    count=count+1;
   
    
    point_ints=line_intersect_plane(xyz(iint,:),xyz(iint+1,:),orient_slice.center,norm);
    
    ijk=xyz2ijk(point_ints,orient_slice);
    
    m(ijk(1),ijk(2),ijk(3))=i;
    
end

fprintf('%d intersecting PVS found',count);

id=unique(m(:));
id(id==0)=[];
disp(id');




