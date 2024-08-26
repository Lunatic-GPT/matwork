
function center=pos2center(pos,vsize,rotmat,sz)

 
    sz=sz(:)';
    pos=pos(:);
    
    ij_c=floor(double(sz(1:2))/2);      
    ijk_c=[ij_c,(sz(3)-1)/2];
    center=pos+(rotmat.*vsize(:)')*ijk_c';
