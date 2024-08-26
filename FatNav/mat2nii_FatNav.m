function mat2nii_FatNav(fmat,fpro,voxsize,center0)
% center0: true: new sequence at VE11; false: sequence at 7 T

pos=zeros(1,3);
fov=zeros(1,3);
norm=zeros(1,3);
if strcmp(fpro(end-3:end),'.pro')
    load(fmat);
    pos(1) = readsPar(fpro,'sCuboid.sPosition.dSag');
    pos(2)= readsPar(fpro,'sCuboid.sPosition.dCor');
    pos(3)=readsPar(fpro,'sCuboid.sPosition.dTra');
    %  fov(1) = readsPar(fullfile(prefix,[prefix,'.pro']),'sCuboid.dReadoutFOV');
    %  fov(2) = readsPar(fullfile(prefix,[prefix,'.pro']),'sCuboid.dPhaseFOV');
    %  fov(3) = readsPar(fullfile(prefix,[prefix,'.pro']),'sCuboid.dThickness');
    if ~exist('voxsize','var')
        o.voxsize=readsPar(pro,'adFree[9]')*ones(1,3);
    else
        o.voxsize=voxsize*[1,1,1];
    end
    o.center=pos;
    
    norm(1)=readsPar(fpro,'sCuboid.sNormal.dSag');
    norm(2)=readsPar(fpro,'sCuboid.sNormal.dCor');
    norm(3)=readsPar(fpro,'sCuboid.sNormal.dTra');
    dphi=readsPar(fpro,'sCuboid.dInPlaneRot');
    if center0
      rotmat=get_rotmat_new(norm,dphi); %new sequence
    else
      rotmat=get_rotmat_old(norm,dphi); %old sequence
    end
    
    o.rotmat=rotmat;
    
    
    o.orient='SPR';
    o.sz=size(d);
    o=center2pos_o(o,o.center);
    
    
    mat2niigz(fmat,'d',strtok(fmat,'.'),false,o);
    
    
else
    a=mapVBVD(fpro);
    load(fmat);
    fov(1)=a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
    fov(2)=a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
    fov(3)=a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
    

   if isfield(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1},'sPosition')
     pos(1)=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition,'dSag');
     pos(2)=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition,'dCor');
     pos(3)=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition,'dTra');   
   else
     pos=[0,0,0];
   end

    if ~exist('voxsize','var')
        o.voxsize=a{2}.hdr.MeasYaps.sWipMemBlock.adFree{10}*ones(1,3);
    else
        o.voxsize=voxsize*[1,1,1];
    end
    o.center=pos;

   dSag=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal,'dSag');
   dCor=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal,'dCor');
   dTra=get_field(a{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal,'dTra');
 
    norm=[dSag,dCor,dTra];
    dphi=0; %don't know how to get the inplane rotation angle
    rotmat=NormInplaneRot2Rotmat(norm,dphi);
    rotmat(:,[2,3])=-rotmat(:,[2,3]);
    o.rotmat=rotmat;  %now correct for AD01
    % assume
    % temporary for Test_Zong; should calculate from header.
    
    %v1=[-0.0366,-0.0018,-0.9993];
    %v2=-[-0.0489,0.9988,-0.0000];
    % v3=cross(v1,v2);
    % o.rotmat=[v1;v2;v3]';
    
    %o.rotmat=[0,0,-1;0,-1,0;-1,0,0]';
    
    
    o.orient='SPR';
    o.sz=size(d);
    o=center2pos_o(o,o.center);
    
    
    mat2niigz(fmat,'d',strtok(fmat,'.'),false,o);
    
end

function rotmat=get_rotmat_new(norm,dphi)

    rotmat=NormInplaneRot2Rotmat(norm,dphi);
    %rotmat(:,2)=-rotmat(:,2); %need to confirm (;,3) left - right;
%tested for PC_MOCO_TEST01; tra; dphi = 0;
  if norm(3)==1 && dphi==0
    rotmat(:,[1,2])=rotmat(:,[2,1]);
  end
    rotmat(:,[2,3])=-rotmat(:,[2,3]);
    
function rotmat=get_rotmat_old(norm,dphi)

%tested for sag, inplanerot = 0
    rotmat=NormInplaneRot2Rotmat(norm,dphi);
    rotmat(:,[2,3])=-rotmat(:,[2,3]);
   
    
function res=get_field(s,fname)
    if isfield(s,fname)
      res=getfield(s,fname);
   else
    res=0;
   end

    
