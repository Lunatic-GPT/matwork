function nii=nii_from_orient(nii,orient,method)
   % It appears ITK-snap does not support method 3 ;
  % method 2: qform_code>0 but sform_code==0; 
  % method 3: sform_code>0;  

if ~exist('method','var')
    method=2;
end
 % orient: should contain voxsize, pos, and rotmat; 
 
 voxsize=orient.voxsize(:)';
 pos=orient.pos;
 rotmat=orient.rotmat;
 rotmat(1:2,:)=-rotmat(1:2,:);
 pos(1:2)=-pos(1:2);
 
   
    nii.untouch=1;
    nii.hdr.hist.magic='n+1';
 % apply both method 2 and method 3    
 %if method==2
    nii.hdr.hist.qform_code=1; %was 0
    nii.hdr.hist.sform_code=0; % was 2

    nii.hdr.hist.qoffset_x=pos(1);
    nii.hdr.hist.qoffset_y=pos(2);
    nii.hdr.hist.qoffset_z=pos(3);
 
   [b,c,d,qfac]=rotmat2quaternion(rotmat);
    nii.hdr.dime.pixdim(1:4)=[qfac,voxsize];
    nii.hdr.hist.quatern_b=b;
    nii.hdr.hist.quatern_c=c;
    nii.hdr.hist.quatern_d=d;
    
 
% else
     nii.hdr.hist.qform_code=1;
    nii.hdr.hist.sform_code=2; 
    nii.hdr.hist.srow_x=[rotmat(1,:).*voxsize,pos(1)];
    nii.hdr.hist.srow_y=[rotmat(2,:).*voxsize,pos(2)];
    nii.hdr.hist.srow_z=[rotmat(3,:).*voxsize,pos(3)];
% end