function write_afni(brik,prefix,info,do_scale)
  % write_afni(brikdata,prefix[,info,do_scale])
  % do_scale: if true data will be scaled and saved as short to save disk
  % space.  Otherwise, saved as float. default: false.
  % the info structure can contain the follow fields
  %  info.ORIGIN. origin position([0,0,0]);  R2L,A2P,I2S
  %  info.DELTA. voxel size (default [1,1,1]);  negative if
  %  different from [0 3 4]
  %  info.ORIENT_SPECIFIC(default [0 3 4]);(R2L,A2P,I2S; afni convention);  
  
  %  info.TAXIS_FLOATS = [0] = Time origin (units: seconds)              This is 0 in datasets created by to3d (at present).
  %             [1] = Time step (TR).
  %             [2] = Duration of acquisition.  This is 0 in datasets
  %                   created by to3d (at present)
  %             [3] = If TAXIS_NUMS[1] > 0, then this is the z-axis offset
  %                   for the slice-dependent time offsets.  This will
  %                   be equal to ORIGIN[2] in datasets created by to3d.c.
  %             [4] = If TAXIS_NUMS[1] > 0, then this is the z-axis step
  %                   for the slice-dependent time offsets.  This will
  %                   be equal to DELTA[2] in datasets created by to3d.c.

  if ~exist('do_scale','var');
      do_scale=false;
  end
  
  if ~exist('info','var')
      info =[];
  end
  
  if ~isfield(info,'ORIGIN')
      info.ORIGIN=[0,0,0];
  end
  if ~isfield(info,'DELTA')
    info.DELTA = [1,1,1];
  end
  if ~isfield(info,'ORIENT_SPECIFIC')
    info.ORIENT_SPECIFIC = [0 3 4];  
  end
  if ~isfield(info,'TAXIS_FLOATS')
    info.TAXIS_FLOATS = [0,2,0,0,0];
  end
  
   info.TYPESTRING = '3DIM_HEAD_ANAT';
   info.SCENE_DATA = [0,2,0];
   info.HISTORY_NOTE = 'write_afni';
   info.DATASET_DIMENSIONS = [size(brik,1),size(brik,2),size(brik,3)];
   info.DATASET_RANK(1) = 3;
   info.DATASET_RANK(2) = size(brik,4);  % the number of subbricks to save.
 
   if do_scale
     brik=round(brik/max(abs(brik(:)))*(2^15-1));
       info.BRICK_TYPES=ones(1,size(brik,4));  %1: short, 3: float
   else
       info.BRICK_TYPES=3*ones(1,size(brik,4));  %1: short, 3: float
   end
       info.BRICK_LABS = [];
       info.IDCODE_STRING = '';
       info.BRICK_STATS=[];  %minimum and maximum values of the subbrick;
       info.BRICK_STATAUX = [];
       info.BRICK_FLOAT_FACS = [];
       info.IDCODE_DATE = date;
       
        info.TAXIS_NUMS=[size(brik,4),0,77002];
         
           %no slice-dependenno slice-dependent time offsets are present (all slices
                   %are presumed to be acquired at the same time)
         % units of time: second
        opt.AdjustHeader='n';
       opt.Prefix = prefix;
       opt.OverWrite = 'y';
       WriteBrik(brik,info,opt);

   