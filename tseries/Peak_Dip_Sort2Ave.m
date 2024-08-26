function Peak_Dip_Sort2Ave(bname,pk_ind,dp_ind,bl,nTR_trial,prefix)
%Peak_Dip_Sort2Ave(bname,pk_ind,dp_ind,bl,nTR_trial,prefix)
% peak and dip positions are specified explicitly.
% pk_ind,dp_ind,bl: 1 based.
% calculate undershoot and peak values in each voxel.
% bname: sorted time course.
% the output data will consist of the following sub-briks:
% 1: peak height
% 2: SEM of peak height
% 3: z value
% 4: peak position
% 5: undershoot depth
% 6: SEM of undershoot depth
% 7: z value
% 8: undershoot position
% 9: ratio between mean undershoot and mean peak height
% 10: mean of ratio between undershoot and peak height 
% 11: SEM of the ratio
% 12: z value

[a,info] = BrikLoadf(bname);

sz = size(a);
if mod(sz(4),nTR_trial) ~=0
    error('mod(sz(4),nTR_trial) ~=0');
end

a_rshp = reshape(a,[sz(1:3),nTR_trial,sz(4)/nTR_trial]);

if ~isempty(bl)
   ts = a_rshp - repmat(mean(a_rshp(:,:,:,bl,:),4),[1,1,1,nTR_trial,1]);
else
    ts = a_rshp;
end

ts_mn = squeeze(mean(ts,5));


data = zeros([sz(1:3),12]);

data(:,:,:,1) = mean(ts_mn(:,:,:,pk_ind),4);
data(:,:,:,5) = mean(ts_mn(:,:,:,dp_ind),4);
data(:,:,:,4) = mean(pk_ind);
data(:,:,:,8) = mean(dp_ind);

pk_arr = squeeze(mean(ts(:,:,:,pk_ind,:),4));
dp_arr = squeeze(mean(ts(:,:,:,dp_ind,:),4));


data(:,:,:,2) = std(pk_arr,0,4)/sqrt(size(pk_arr,4));
t1 = data(:,:,:,1)./data(:,:,:,2);
data(:,:,:,6) = std(dp_arr,0,4)/sqrt(size(dp_arr,4));
t2 = data(:,:,:,5)./data(:,:,:,6);
t1(isnan(t1)) = 0;
t2(isnan(t2)) = 0;

data(:,:,:,9) = data(:,:,:,5)./data(:,:,:,1);
data(:,:,:,10) = mean(dp_arr./pk_arr,4);
data(:,:,:,11) = std(dp_arr./pk_arr,0,4)/sqrt(size(pk_arr,4));
t3 = data(:,:,:,10)./data(:,:,:,11);
t3(isnan(t3)) = 0;

data(isnan(data)) =  0;
for i=1:sz(1)
    for j=1:sz(2)
        for k=1:sz(3)      
            if sz(4)/nTR_trial>1
               data(i,j,k,3) = t2z(t1(i,j,k),sz(4)/nTR_trial-1);
               data(i,j,k,7) = t2z(t2(i,j,k),sz(4)/nTR_trial-1);
               data(i,j,k,12) = t2z(t3(i,j,k),sz(4)/nTR_trial-1);
            end
        end
    end
end

if exist('prefix','var')
       
       info.IDCODE_STRING = 'Peak_Dip';
       
       info.TYPESTRING = '3DIM_HEAD_ANAT';
       info.SCENE_DATA = [0 2 0];  % first 0: orig view; second 2: ANAT_EPI_TYPE; 0: matches TYPESTRING 
       
       info.BRICK_STATS=[];  %minimum and maximum values of the subbrick;
       
       
       info.BRICK_LABS=['Peak Height~SEM of Height~z score~Peak Position~'...
                        'Undershoot depth~SEM of undershoot~z score~Dip Position~'...
                        'ratio of mean us to mean pk~mean of ratios of us to pk~' ...
                        'SEM of ratio~z score~'];
       
       info.BRICK_FLOAT_FACS = [];
       
       info.IDCODE_DATE = date;
       info.TAXIS_NUMS=[];
       
       history1 = sprintf('Peak_Dip_Sort2Ave(%s,%s,%s,%s,%d,%s)',...
           bname,num2str(pk_ind),num2str(dp_ind),num2str(bl),nTR_trial,prefix);
     
       info.HISTORY_NOTE = [info.HISTORY_NOTE,'\n',history1];

       nSb = size(data,4);
       info.DATASET_RANK(2) = nSb;  % the number of subbricks to save.
       info.BRICK_TYPES= 3*ones(1,nSb);  %1: short, 3: float
             
       opt.Prefix = prefix;
       opt.OverWrite = 'y';
       
       WriteBrik(data,info,opt);
       
end
