function do_QSM_star(dcm_ph,steps,skull_strip,outpath,thk)
%do_QSM_star(dcm_ph,steps[,skull_strip,outpath])
% steps: 1 only do phase unwrapping 
% steps: 2  background phase removal
% steps: 4 do QSM

% thickness in mm
%%
if ~exist('skull_strip','var')

    skull_strip='skull_strip_out+orig.HEAD';
end
mask= BrikLoad(skull_strip);  %% load mask
mask=mask>0;

mask=flip(mask,3);

TE = readdPar(dcm_ph,'EchoTime');% load TE in ms



a=readdPar(dcm_ph,'ImageOrientationPatient');
spacing=readdPar(dcm_ph,'PixelSpacing');
if ~exist('thk','var')
thk=readdPar(dcm_ph,'SliceThickness');
end

b = cross(a(1:3), a(4:6));
H = [ a(6) a(3) b(3)];


voxelsize = [spacing' thk];


ph=ri(dcm_ph,[]);
%mag=ri(dcm_mag,[]);

ph=permute(ph,[2,1,3]);

ph=ph*2*pi/4096-pi;


padsize = [12 12 12];


  %Unwrapped_Phase=LaplacianPhaseUnwrap(ph,'voxelsize',voxelsize,'padsize',padsize);
  Unwrapped_Phase=MRPhaseUnwrap(ph,'voxelsize',voxelsize,'padsize',padsize);
%   
  if ~exist('outpath','var') || isempty(outpath)
      outpath='';
      cur_dir=pwd;
  else
      mkdir(outpath);
      cur_dir=cd(outpath);
  end
  
save_nii(make_nii(Unwrapped_Phase,voxsize),'UwPhase.nii.gz');


  [TissuePhase,UpdatedMask] =V_SHARP(single(Unwrapped_Phase),single(mask),'smvsize',25,'voxelsize',voxelsize);
save_nii(make_nii(UpdatedMask),'UpdatedMask.nii.gz');
save_nii(make_nii(TissuePhase),'TissuePhase.nii.gz');

if steps==2
    star = QSM_star(single(TissuePhase) ,single(UpdatedMask),'H',H,'voxelsize',voxelsize,'padsize',padsize,'TE',TE,'B0',3);  
    save_nii(make_nii(star,voxelsize),'QSM.nii.gz');
end

cd(cur_dir);

