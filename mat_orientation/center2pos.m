function pos=center2pos(vsize,rotmat,sz,center)
% pos=center2pos(vsize,rotmat,sz,center)
if length(sz)==2
    sz(3)=1;
end

    ij_c=floor(double(sz(1:2))/2);      
    ijk_c=[ij_c,(sz(3)-1)/2];
    pos=center(:)-(rotmat.*vsize(:)')*ijk_c';

    
    
   