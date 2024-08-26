function Peak_Dip_Sort2Ave3(bname,nTR_pk, nTR_dp,bl,nTR_trial,prefix)
%Peak_Dip_Sort2Ave3(bname,nTR_pk,nTR_dip,bl,nTR_trial,prefix)
% peak positions are determined automatically by finding the
% maximum by summing over the neighboring nTR_pk time points.
% undershoot is assumed to start from the first negative data point post-peak.
% bname: sorted time course.
% bl: 1 based.
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

[a,info] = BrikLoad(bname);

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


ind_nbl = setdiff(1:nTR_trial,bl);
ts_pk = zeros([sz(1:3),length(ind_nbl)-nTR_pk+1]);


for i=1:length(ind_nbl)-nTR_pk+1
  i1 = ind_nbl(i);  
  ts_pk(:,:,:,i) = mean(ts_mn(:,:,:,i1:i1+nTR_pk-1),4);
end


data = zeros([sz(1:3),12]);

[data(:,:,:,1),ipk] = max(ts_pk,[],4);

idp = zeros(sz(1:3));
for i=1:sz(1)
    for j=1:sz(2)
        for k=1:sz(3)
          
            tmp = find(ts_mn(i,j,k,ipk(i,j,k):end)<=0,1);
            if isempty(tmp)
                continue;
            end
            
            idp(i,j,k) = tmp+ipk(i,j,k)-1;
            idp_tmp = idp(i,j,k);
            if isempty(intersect(bl,idp_tmp:idp_tmp+nTR_dp-1)) && idp_tmp+nTR_dp-1<=nTR_trial
             data(i,j,k,5) = mean(ts_mn(i,j,k,idp_tmp:idp_tmp+nTR_dp-1),4);
            end
        end
    end
end


data(:,:,:,4) = ipk+(nTR_pk-1)/2;
data(:,:,:,8) = idp+(nTR_dp-1)/2;



pk_arr = zeros([sz(1:3),sz(4)/nTR_trial]);
dp_arr = zeros([sz(1:3),sz(4)/nTR_trial]);

for i=1:sz(1)
    for j=1:sz(2)
        for k = 1:sz(3)
            
          pi = ipk(i,j,k);
          di = idp(i,j,k);
          pk_arr(i,j,k,:) = mean(ts(i,j,k,pi:pi+nTR_pk-1,:),4);
          if di >0 && di+nTR_dp-1<=nTR_trial && isempty(intersect(bl,di:di+nTR_dp-1))
            dp_arr(i,j,k,:) = mean(ts(i,j,k,di:di+nTR_dp-1,:),4);
          end
        end
    end
end

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


for i=1:sz(1)
    for j=1:sz(2)
        for k=1:sz(3)      
            if sz(4)/nTR_trial >1
              data(i,j,k,3) = t2z(t1(i,j,k),sz(4)/nTR_trial-1);
              data(i,j,k,7) = t2z(t2(i,j,k),sz(4)/nTR_trial-1);
              data(i,j,k,12) = t2z(t3(i,j,k),sz(4)/nTR_trial-1);
            end
        end
    end
end


data(isnan(data)) =  0;

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
       history1 = sprintf('Peak_Dip_Sort2Ave2(%s,%d,%d,%s,%d,%s)',...
           bname,nTR_pk,nTR_dp,num2str(bl),nTR_trial,prefix);
     
       info.HISTORY_NOTE = [info.HISTORY_NOTE,'\n',history1];

       nSb = size(data,4);
       info.DATASET_RANK(2) = nSb;  % the number of subbricks to save.
       info.BRICK_TYPES= 3*ones(1,nSb);  %1: short, 3: float
             
       opt.Prefix = prefix;
       opt.OverWrite = 'y';
       
       WriteBrik(data,info,opt);
       
end