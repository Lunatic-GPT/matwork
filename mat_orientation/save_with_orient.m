function save_with_orient(fname,orient,d)


voxsize=orient.voxsize;
center=orient.center;
rotmat=orient.rotmat;
pos=orient.pos;
orient=orient.orient;

save(fname,'d','voxsize','center','rotmat','pos','orient');
