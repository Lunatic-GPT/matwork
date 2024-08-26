function ijk=xyz2ijk(xyz,orient)
% pos=center2pos(vsize,rotmat,sz,center)

  pos=orient.pos;
  rotmat=orient.rotmat;
  vsize=orient.voxsize;
  
  ijk=(rotmat.*vsize(:)')\(xyz(:)-pos(:));  
  ijk=round(ijk)+1;
  
  
  

    
   