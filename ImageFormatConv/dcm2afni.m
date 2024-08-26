function dcm2afni(dname,pattern,nslices,prefix)
tic;

     str = dir(fullfile(dname,pattern));

     if mod(length(str),nslices) ~=0
         error('number of slices parameter wrong');
     end
     
     nt = length(str)/nslices;
     
    for i=1:length(str)
        
        img = dicomread(fullfile(dname,str(i).name));
        if i==1
            dinfo = dicominfo(fullfile(dname,str(i).name));
            brikdata = zeros([size(img),nslices,nt]);
        end
        
        brikdata(:,:,mod(i-1,nslices)+1,ceil(i/nslices)) = img;
    end
    
        


       info.DATASET_RANK = [3,nt,0,0,0,0,0,0];  % the number of subbricks to save.
       info.BRICK_TYPES= 3*ones(1,nt);  %1: short, 3: float

       sz = size(brikdata);
       info.DATASET_DIMENSIONS = [sz(1:3),0,0];
       info.BRICK_LABS = [];
       info.IDCODE_STRING = 'dcm2afni images';
       info.ORIENT_SPECIFIC = [0,3,4];
       info.ORIGIN = [0 0 0];
       info.DELTA = [1.7188 1.7188 3];
       info.BYTEORDER_STRING='LSB_FIRST';
       info.TYPESTRING = '3DIM_HEAD_ANAT';
       info.SCENE_DATA = [0 2 0];  % first 0: orig view; second 2: ANAT_EPI_TYPE; 0: matches TYPESTRING 
       
       info.BRICK_STATS=[];  %minimum and maximum values of the subbrick;
       
       info.BRICK_FLOAT_FACS = [];
       
       info.IDCODE_DATE = date;
       info.TAXIS_NUMS=[nt nslices 77002 -999 -999 -999 -999 -999];
       info.TAXIS_FLOATS=[0 90 0 0 0 -999999 -999999 -999999];
       info.TAXIS_OFFSETS = zeros(1,nslices);
       info.HISTORY_NOTE = sprintf('dcm2afni(%s,%s,%d)',dname,pattern,nslices);
       opt.OverWrite = 'y';
       opt.Prefix = prefix;
       
       WriteBrik(brikdata,info,opt);
         
       
     
%% save the average time series for each trial type.    

    disp([mfilename ' finish in ', num2str(toc), ' s']);