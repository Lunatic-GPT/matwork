function xyz=ijk2xyz(ijk,orient)
% pos=center2pos(vsize,rotmat,sz,center)
% ijk:n*3
% xyz: n*3;
  pos=orient.pos;
  rotmat=orient.rotmat;
  vsize=orient.voxsize;
  
  xyz=pos(:)+(rotmat.*vsize(:)')*(ijk'-1);
  xyz=xyz';
  
  
  
  
  
  

    
   