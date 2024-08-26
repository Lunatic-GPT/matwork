function [an,ipvs]=find_PVS_intersect_slice(pvs,orient,center,norm)
% input:
% pvs: the structure of PVSs
% orient: the orientation of PVS
% center and norm: the center and normal direction of the PC-MRI slice
% an: angles of the PVS that intersect the plane
% ipvs: indices of the associated PVS;


an=[];
ipvs=[]; % the indices of the PVS that intersect the plane.
count=0;

for i=1:length(pvs)
    sub=pvs(i).subind(1:pvs(i).nvox_path,:);
    
    xyz=ijk2xyz(sub,orient);
    res=dist2plane(xyz,center,norm);
    
    
    iint=find(sign(res(1:end-1)).*sign(res(2:end))<=0);
    if length(iint)>1 || isempty(iint)
        continue;
    end
    
    count=count+1;
    if abs(res(iint))<abs(res(iint+1))
        iint2=iint;
    else
        iint2=iint+1;
    end
    an(count)=angle_bw_2vec(pvs(i).dir_seg(:,iint2)',norm);
    ipvs(count)=i;
      
end
