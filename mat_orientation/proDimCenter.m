function coord=proDimCenter(protocol,sz)
%proDimCenter(dname) or dcmDimCenter(pos,ImageOrientPatient,voxsize)
% only tested for 2D recon_fl_fq images;
% pos: 3*nslice; 
% orient: 1*6; ImageOrientationPatient; columns then rows, then slice
% voxsize: 1*2; pixsize from dicom header; for column (2nd dim) and then row
%

% output:
% vsize: 1*3; rows then columns; same for orient label 
% center: 1*3
% orient: label
% pos: 3*1; position of pixel 1,1,1
% rotmat: 3*3 matrix; physical position is rotmat*([i-1,j-1,k-1].*voxsize)'+pos; dicom convention;
% % center of the image: for the first two dimensions the center is on the
% n/2+1 (even n) or (n+1)/2 (odd n) pixels, and for third dimension on
% (n+1)/2, regardless of n;  This is the same as the Siemens scanner
% convention and true for both 2D and 3D sequences.
% dicom convention: 
%  x=R-L, y=A-P, and z=I-S (R,A,I are < 0; L,P,S are > 0).
% nifti convention: L-R; P-A; I-S; 


%% this function assume single oblique for now.

% center=get_center(pos,voxsize,OrientPat,columns,rows);
% 
% orient=orientLabel(reshape(vori,[3,2])',pos);
%  


%%
rotmat=NormInplaneRot2Rotmat(protocol); 

rotmat(:,[1,2])=rotmat(:,[2,1]); 
rotmat(:,2)=-rotmat(:,2);
warning('rotmat only tested for PE along PA');

center(1)=readsPar(protocol,'asSlice[0].sPosition.dSag');
center(2)=readsPar(protocol,'asSlice[0].sPosition.dCor');
center(3)=readsPar(protocol,'asSlice[0].sPosition.dTra');

voxsize(3)=readsPar(protocol,'asSlice[0].dThickness');
fov(1)=readsPar(protocol,'asSlice[0].dReadoutFOV');
fov(2)=readsPar(protocol,'asSlice[0].dPhaseFOV');
fov(2)=readsPar(protocol,'asSlice[0].dPhaseFOV');

nvox(1)=readsPar(protocol,'lBaseResolution');
nvox(2)=readsPar(protocol,'lPhaseEncodingLines');

if ~exist('sz','var')
    sz=[nvox,1];
end

if length(sz)==2
    sz(3)=1;
end
voxsize(1:2)=fov./sz(1:2);
orient=orientLabel(rotmat);

pos=center2pos(voxsize,rotmat,sz,center);

coord.voxsize=voxsize;
coord.center=center;
coord.rotmat=rotmat;
coord.pos=pos;
coord.orient=orient;


function res=get_rotmat(protocol)
%%
InplaneRot=readsPar(protocol,'asSlice[0].dInPlaneRot');
% 
% [th,ph]=unitVec2thetaPhi(norm');  %
% rotmat_norm=transform_matrix_rotation(th,ph+90);

norm(1)=readsPar(protocol,'asSlice[0].sNormal.dSag');
norm(2)=readsPar(protocol,'asSlice[0].sNormal.dCor');
norm(3)=readsPar(protocol,'asSlice[0].sNormal.dTra');

norm_nz=norm(abs(norm)>1e-6);
if length(norm_nz)>2
    error('double oblique not implemented');
end


if abs(norm(2))<1e-6 % S-T
    rotmat_norm=[norm(3),0,norm(1);0,-1,0;-norm(1),0,norm(3)];
rotmat_inplane=[cos(InplaneRot),-sin(InplaneRot),0;sin(InplaneRot),cos(InplaneRot),0;0,0,1];
res=rotmat_norm*rotmat_inplane;

elseif abs(norm(1))<1e-6 % S-C
   
   % rotmat_norm=[-1,0,0;0,norm(3),norm(2);0,-norm(2),norm(3)];
   rotmat_norm=[1,0,0;0,-norm(3),norm(2);0,norm(2),norm(3)];
   
    rotmat_inplane=[cos(InplaneRot),-sin(InplaneRot),0;sin(InplaneRot),cos(InplaneRot),0;0,0,1];
    res=rotmat_norm*rotmat_inplane;

end


% 
% a=setdiff(1:3,f2);
% if length(pos)>1
%     dpos=pos(2)-pos(1);
% else
%     dpos=1;  %does not matter
% end
% 
% j=a*2+(sign(dpos)-1)/2;
% ori(3)=label(j);
% 

  
  



