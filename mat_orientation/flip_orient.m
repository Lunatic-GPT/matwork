function res=flip_orient(o,dim)

for i=1:length(dim)

    o.rotmat(:,dim(i))=-o.rotmat(:,dim(i));
    o.pos=o.pos-o.rotmat(:,dim(i))*(o.sz(dim(i))-1)*o.voxsize(dim(i));

end

o.center=pos2center(o.pos,o.voxsize,o.rotmat,o.sz);
res=o;