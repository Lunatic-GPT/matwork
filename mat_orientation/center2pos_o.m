function o=center2pos_o(o,center)
% pos=center2pos(vsize,rotmat,sz,center)
sz=o.sz;
vsize=o.voxsize;
rotmat=o.rotmat;

if length(sz)==2
    sz(3)=1;
end

    ij_c=floor(double(sz(1:2))/2);      
    ijk_c=[ij_c,(sz(3)-1)/2];
    o.pos=center(:)-(rotmat.*vsize(:)')*ijk_c';
    o.center=center;
    
    
    
   