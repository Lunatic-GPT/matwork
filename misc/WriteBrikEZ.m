function WriteBrikEZ(brik,info,history,prefix,labels)
% WriteBrikEZ(brik,info,history,prefix,labels)
       
if ~exist('labels','var')
    labels = [];
end

       info.DATASET_DIMENSIONS = [size(brik,1),size(brik,2),size(brik,3)];
       info.DATASET_RANK = [3, size(brik,4)];  % the number of subbricks to save.
       info.BRICK_TYPES=3*ones(1,size(brik,4));  %1: short, 3: float
       info.BRICK_LABS = labels;
       info.IDCODE_STRING = '';
       info.BRICK_STATS=[];  %minimum and maximum values of the subbrick;
       info.BRICK_STATAUX = [];
       info.BRICK_FLOAT_FACS = [];
       info.IDCODE_DATE = date;
       info.TYPESTRING='3DIM_HEAD_ANAT';  %it should match scene_data(3)  
         info.SCENE_DATA(3)=0;
         
         info.TAXIS_NUMS(1)=size(brik,4);
         info.TAXIS_NUMS(2)=0;  %no slice-dependenno slice-dependent time offsets are present (all slices
                   %are presumed to be acquired at the same time)
         info.TAXIS_NUMS(3)=77002;  % units of time: second
         
         info.TAXIS_FLOATS(1) = 0;
         info.TAXIS_FLOATS(2) = 1;
         info.TAXIS_FLOATS(3) = 0;
         info.TAXIS_FLOATS(4) = 0;
         info.TAXIS_FLOATS(5) = 0;
         if isfield(info,'HISTORY_NOTE')
         info.HISTORY_NOTE = [info.HISTORY_NOTE, '\n',history];
         else
          info.HISTORY_NOTE = history;
         end
       opt.AdjustHeader='n';
       opt.Prefix = prefix;
       opt.OverWrite = 'y';
       WriteBrik(brik,info,opt);
