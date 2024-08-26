function write_afni_sdt(brik,info,prefix,format)
  % write_afni(brikdata,delta,prefix)
  % the info structure can contain the follow fields
  %  info.ORIGIN. origin position([0,0,0]);
  %  info.DELTA. voxel size (default [1,1,1]);
  %  info.ORIENT_SPECIFIC(default [1 3 5]);(L2R,A2P,S2I)  (left to right,
  %  up to down, in to out).
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

  if strcmp(format,'a')
      
  if ~isfield(info,'ORIGIN')
      info.ORIGIN=[0,0,0];
  end
  if ~isfield(info,'DELTA')
    info.DELTA = [1,1,1];
  end
  if ~isfield(info,'ORIENT_SPECIFIC')
    info.ORIENT_SPECIFIC = [1 3 5];
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
   info.BRICK_TYPES=3*ones(1,size(brik,4));  %1: short, 3: float
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
       opt.Prefix = prefix;
       opt.OverWrite = 'y';
       WriteBrik(brik,info,opt);

  elseif strcmp(format,'s')
   writesdt4(brik,prefix);    
  else
      error('unknown format');
  end
  