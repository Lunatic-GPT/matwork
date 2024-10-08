function Peak_Dip(bname,bl,nTR_pk,nTR_dp,method,prefix)
%[us,pk] = Peak_Dip(bname,bl[,prefix])
% bl: 1 based.
% calculate undershoot and peak values in each voxel.
% method: the method for calculating undershoot. 
%        1: minimum of the mean of neighbouring nTR_us points.
%        2: Mean of nTR_us time points following positive response.  


[a,info] = BrikLoad(bname);


a = a - repmat(mean(a(:,:,:,bl),4),[1,1,1,size(a,4)]);
sz = size(a);

ind_nbl = setdiff(1:sz(4),bl);

ts_pk = zeros([sz(1:3),length(ind_nbl)-nTR_pk+1]);

d = zeros([sz(1:3),2]);

for i=1:length(ind_nbl)-nTR_pk+1
  i1 = ind_nbl(i);  
  ts_pk(:,:,:,i) = mean(a(:,:,:,i1:i1+nTR_pk-1),4);
end

[d(:,:,:,1),i_pk] = max(ts_pk,[],4);

if method == 1
    
 ts_dp = zeros([sz(1:3),length(ind_nbl)-nTR_dp+1]);
 for i=1:length(ind_nbl)-nTR_dp+1
   i1 = ind_nbl(i);
   ts_dp(:,:,:,i) = mean(a(:,:,:,i1:i1+nTR_dp-1),4);
 end

else
   
    for i=1:sz(1)
        for j=1:sz(2)
            for k=1:sz(3)
                i_start = ind_nbl(i_pk(i,j,k));
                i_us_first = find(a(i,j,k,i_start:ind_nbl(end))<0,1);
                
                if ~isempty(i_us_first)
                    i_us_first = i_us_first + i_start-1;
                    d(i,j,k,2) = mean(a(i,j,k,i_us_first:i_us_first+nTR_dp-1));
                end
            end
        end
    end
  
    
end   

if exist('prefix','var')
       
       info.IDCODE_STRING = 'Peak_Dip';
       info.TYPESTRING = '3DIM_HEAD_ANAT';
       info.SCENE_DATA = [0 2 0];  % first 0: orig view; second 2: ANAT_EPI_TYPE; 0: matches TYPESTRING 
       
       info.BRICK_STATS=[];  %minimum and maximum values of the subbrick;
       info.BRICK_LABS='peak~undershoot~';  %minimum and maximum values of the subbrick;
       
       info.BRICK_FLOAT_FACS = [];
       
       info.IDCODE_DATE = date;
       info.TAXIS_NUMS=[];
       history1 = sprintf('Peak_Dip(%s,%s,%d,%d,%d,%s)',bname,num2str(bl),nTR_pk,nTR_dp,method,prefix);
       info.HISTORY_NOTE = [info.HISTORY_NOTE,'\n',history1];

      
       nSb = size(d,4);
       info.DATASET_RANK(2) = nSb;  % the number of subbricks to save.
       info.BRICK_TYPES= 3*ones(1,nSb);  %1: short, 3: float
             
       opt.Prefix = prefix;
       opt.OverWrite = 'y';
       
       WriteBrik(d,info,opt);
       
end
