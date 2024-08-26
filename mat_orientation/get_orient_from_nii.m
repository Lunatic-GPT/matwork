
function orient=get_orient_from_nii(nii,dim_flip)
% dicom convention; which has x and y axes reversed compared with nifti
% convention;
% dim_flip: the dimension(s) to flip
if ~exist('dim_flip','var')
    dim_flip=[];
end

if isa(nii,'char')
    nii=load_untouch_niigz(nii);
end

voxsize=nii.hdr.dime.pixdim(2:4);

if nii.hdr.hist.sform_code==0 && nii.hdr.hist.qform_code>0% method 2
    rotmat=zeros(3,3);
    qfac=nii.hdr.dime.pixdim(1);
    if qfac==0
        qfac=1;
    end
    
    b=nii.hdr.hist.quatern_b;
    c=nii.hdr.hist.quatern_c;
    d=nii.hdr.hist.quatern_d;
    a=sqrt(1-b*b-c*c-d*d);
    
    rotmat(1,:)=-[ a*a+b*b-c*c-d*d,  2*b*c-2*a*d,       2*b*d+2*a*c];
    rotmat(2,:)=-[ 2*b*c+2*a*d,       a*a+c*c-b*b-d*d,   2*c*d-2*a*b];
    rotmat(3,:)=[ 2*b*d-2*a*c,       2*c*d+2*a*b,       a*a+d*d-c*c-b*b];
    
    rotmat(:,3)=rotmat(:,3)*qfac;
    pos=[-nii.hdr.hist.qoffset_x;-nii.hdr.hist.qoffset_y;nii.hdr.hist.qoffset_z];
    
else  %method 3
    rotmat=[-nii.hdr.hist.srow_x(1:3);-nii.hdr.hist.srow_y(1:3);nii.hdr.hist.srow_z(1:3)]./repmat(voxsize,[3,1]);
    
    
    pos=[-sum(nii.hdr.hist.srow_x(4));-sum(nii.hdr.hist.srow_y(4));sum(nii.hdr.hist.srow_z(4))]; %pos is the position of the 
end



orient.voxsize=nii.hdr.dime.pixdim(2:4);% not implemented yet

orient.rotmat = rotmat;
orient.pos=pos;
orient.orient=orientLabel(orient.rotmat);

orient.sz=nii.hdr.dime.dim(2:4);

orient.center=pos2center(orient.pos,orient.voxsize,orient.rotmat,orient.sz);


orient=flip_orient(orient,dim_flip);


