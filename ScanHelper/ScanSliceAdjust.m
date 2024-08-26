function res=ScanSliceAdjust(ref,nav1,fpos0,fpos1,for_rev,dfile)
% res=ScanSliceAdjust(ref,nav1,fpos0,fpos1[,for_rev,dfile])
% for_inv: forward (fpos0 is the slice location on ref) or reverse (fpos0 is
% the slice location of a fl_fq scan).
% dfile: default - sprintf('dfile_%s.1D',nav1);

if ~exist('for_rev','var')
    for_rev='forward';
end

sadj=tic;

if ~ispc
    tic;
    if ~exist([ref,'.nii.gz'],'file')
        if ~isempty(dir([ref,'/*.dcm']))
            cmd=sprintf('to3d -prefix %s.nii.gz %s/*.dcm',ref,ref);
        else
            cmd=sprintf('to3d -prefix %s.nii.gz %s/*.IMA',ref,ref);
        end
        system(cmd);  %about 2 s
        toc;
    end
    
    if ~exist([nav1,'.nii.gz'],'file')
        if ~isempty(dir([nav1,'/*.dcm']))
            cmd=sprintf('to3d -prefix %s.nii.gz %s/*.dcm',nav1,nav1);
        else
            cmd=sprintf('to3d -prefix %s.nii.gz %s/*.dcm',nav1,nav1);
        end
        system(cmd);  %about 2 s
        toc;
    end
    if ~exist(sprintf('%s_rs.nii.gz',ref),'file')
        cmd=sprintf('3dresample -dxyz 2 2 2 -prefix %s_rs.nii.gz -input %s.nii.gz',ref,ref);
        system(cmd);  %3 s
    end
    
    if  ~exist([nav1,'_rs.nii.gz'],'file')
        cmd=sprintf('3dresample -master %s_rs.nii.gz -prefix %s_rs.nii.gz -input %s.nii.gz',ref,nav1,nav1);
        system(cmd);
        toc;
    end
    
    cmd=sprintf('3dvolreg -prefix %s_rs_vr.nii.gz -dfile dfile_%s.1D -base %s_rs.nii.gz %s_rs.nii.gz',nav1,nav1,ref,nav1);
    system(cmd);
    
    % tic;
    %      cmd=sprintf('3dresample -master %s.nii.gz -prefix %s_rs.nii.gz -input %s.nii.gz',nav1,ref,ref);
    %      system(cmd);
    %      toc;
    %      tic;
    %      cmd=sprintf('3dvolreg -prefix %s_rs_vr.nii.gz -dfile dfile_%s.1D -base %s_rs.nii.gz %s.nii.gz',nav1,nav1,ref,nav1);
    %     system(cmd);
    %
    %     toc;
end

if ~exist('dfile','var')
    dfile=sprintf('dfile_%s.1D',nav1);
end

res=load(dfile);
m=afni_motionPar2Mat(dfile);
m=reshape(m,[4,3])';
m(4,:)=[0,0,0,1];
mi=inv(m);
mi=mi(1:3,:);
%%
fpos0=str2cell(fpos0);
fpos1=str2cell(fpos1);


for i=1:length(fpos0)
    
    if ~exist(fpos0{i},'file')
        continue;
    end
    
    [pos,norm,inplaneRot]=get_slicePosition(fpos0{i});
    pos=pos(:)';
    
    
    mat=NormInplaneRot2Rotmat(norm,inplaneRot);
    % ctse=FOVcenter_dcm2pro(num2str(ref));
    ctse=get_orient_from_nii([ref,'_rs.nii.gz']);
    ctse=ctse.center(:)';
    %  ctse=0;
    pos=pos-ctse;
    
    if strcmp(for_rev,'forward')
        mat2=mi(:,1:3)*mat;
        pos2=mi*[pos,1]'+ctse';
    else
        mat2=m(1:3,1:3)*mat;
        pos2=m(1:3,:)*[pos,1]'+ctse';
    end
    if exist(fpos1{i},'file')
        movefile(fpos1{i},unique_name(fpos1{i}));
    end
    
    
    rad=Rotmat2InplaneRot(mat2,inplaneRot);
    
    save_slicePosition(fpos1{i},pos2,mat2(:,3),rad);
    
end

toc(sadj);
function res=unique_name(fname)

[prefix,suf]=strtok2(fname,'.');

count='a'-1;
prefix2=prefix;

while 1
    if exist([prefix2,suf],'file')
        count=count+1;
        prefix2 =[prefix,'_',count];
    else
        break;
    end
    
    
end
res=[prefix2,suf];

function [pre,suf]=strtok2(str,tok)
% backward strtok

str=str(end:-1:1);
[suf,pre]=strtok(str,tok);

if ~isempty(pre)
    pre=pre(end:-1:2);
    suf=[tok,suf(end:-1:1)];
else
    [pre,suf]=strtok(str(end:-1:1),tok);
end


function str = str2cell(str)

if ~iscell(str)
    temp = str;
    str = cell(1,1);
    str{1} = temp;
end



function m_all=afni_motionPar2Mat(dfile)
% n  roll  pitch  yaw  dS  dL  dP  rmsold rmsnew
%                 roll  = rotation about the I-S axis }
%                    pitch = rotation about the R-L axis } degrees CCW
%                    yaw   = rotation about the A-P axis }
%                      dS  = displacement in the Superior direction  }
%                      dL  = displacement in the Left direction      } mm
%                      dP  = displacement in the Posterior direction }
% so the motion parameters are the amount of motion from new to base.
% dicom convention; right-hand coordinate system
% the x-axis is increasing to the left hand side of the patient.
% The y-axis is increasing to the posterior side of the patient.
% The z-axis is increasing toward the head of the patient.
% #define ORI_R2L_TYPE  0  /* Right to Left         */
% #define ORI_L2R_TYPE  1  /* Left to Right         */
% #define ORI_P2A_TYPE  2  /* Posterior to Anterior */
% #define ORI_A2P_TYPE  3  /* Anterior to Posterior */
% #define ORI_I2S_TYPE  4  /* Inferior to Superior  */
% #define ORI_S2I_TYPE  5  /* Superior to Inferior  */
% the AFNI convention is that R-L, A-P, and I-S are
%         negative-to-positive. same as DICOM.
% RL 1 dim: AP: 2nd dim; IS: 3rd dim
% the output should be reshaped to ([4,3]);
% the output is matrix (m) is defined such that x=m*x', where x and x' are the original and
% new coordinates of the same point of the object.
if isa(dfile,'char')
    d=load(dfile);
else
    d=dfile;
end

d=d(:,2:7);
m_all=zeros(size(d,1),12);
for i=1:size(d,1)
    
    m_IS = transform_matrix_rotation_arb_axis(0,0,d(i,1));
    m_RL = transform_matrix_rotation_arb_axis(90,0,d(i,2));
    m_AP = transform_matrix_rotation_arb_axis(90,90,d(i,3));
    
    %m=m_IS*m_RL*m_AP;
    %m=m_IS*m_AP*m_RL;
    %m=m_AP*m_RL*m_IS;
    %m=m_AP*m_IS*m_RL;
    %m=m_RL*m_IS*m_AP;
    m=m_AP*m_RL*m_IS;  % not sure which one is correct
    
    % shift=m*d(i,[5,6,4])';
    m(:,end+1)=d(i,[5,6,4]);
    
    m_all(i,:)=vec(m');
end
function a=vec(a)
a=a(:);

function res=transform_matrix_rotation_arb_axis(theta,phi,deg)
% rotating by deg degrees around an axis defined by the angles theta and phi

% the following is ok too
% zp=thetaPhi2unitVec(theta,phi);
%
% [tmp,ind]=max(zp);
% yp=ones(1,3);
% yp(ind)=0;
%
% yp=yp/sqrt(2);
% yp=cross(yp,zp);
% yp=yp/sos(yp,2);
% xp=cross(yp,zp);
%
% m1=[xp',yp',zp'];

m1=transform_matrix_rotation(theta,phi);
m2=transform_matrix_rotation(0,deg);

res=m1*m2/m1;


function res=transform_matrix_rotation(theta,phi)
% theta goes from positive z to positive x
% phi goes from positive x to positive y
% right-handed coordinate system
% angles in degrees
% rotate coordinates first by theta around y; then by phi around z;
% matrix should be on the left, i.e. vec_new=A*vec;

theta=theta*pi/180;
phi=phi*pi/180;

m1=eye(3);

m1([1,3],[1,3])=[cos(theta),sin(theta);-sin(theta),cos(theta)];

m2=eye(3);

m2(1:2,1:2)=[cos(phi),-sin(phi);sin(phi),cos(phi)];

res=m2*m1;

function res3=readsPar(fname,par)
%res=readbPar(fname,par,isnum)

fid=fopen(fname,'r');
res3=[];
while ~feof(fid)
    
    b=fgetl(fid);
    %{
  if b==-1
      fclose(fid);
      if isempty(res3)
        error([par, ' not found']);
      
      else
       return;
      end
  end
    %}
    
    ind=strfind(b,par);
    if ~isempty(ind)
        res=b(ind+length(par):end);
        res=strrm(res,' ');
        res=strrm(res,'=');
        
    else
        continue;
    end
    
    
    res2=str2double(res);
    
    if ~isnan(res2)
        res3(end+1)=res2;
    else
        res3{end+1}=res;
    end
end
if isempty(res3) && nargout==0
    disp([par, ' not found']);
end

fclose(fid);

if isempty(res3)
    res3=0;
end
