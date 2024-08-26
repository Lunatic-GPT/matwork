function nii=nii_new_center_rotmat(nii,o_new)
% assume method = 3 for nifit header coordinates
% get new o_new from 

voxsize=nii.hdr.dime.pixdim(2:4);

nii.hdr.hist.srow_x(1:3)=-o_new.rotmat(1,:).*voxsize;
nii.hdr.hist.srow_y(1:3)=-o_new.rotmat(2,:).*voxsize;
nii.hdr.hist.srow_z(1:3)=o_new.rotmat(3,:).*voxsize;


o=get_orient_from_nii(nii);

o=center2pos_o(o,o_new.center);


nii.hdr.hist.srow_x(4)=-o.pos(1);
nii.hdr.hist.srow_y(4)=-o.pos(2);
nii.hdr.hist.srow_z(4)=o.pos(3);


