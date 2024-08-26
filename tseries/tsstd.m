function tsstd(fnames,tps,out_name)
% tsstd(fnames,tps,out_name)
% fnames: dataset names
% tps: time points included in standard deviation calculation. 1 based. use
% [] to include all time points.
% same for all files.
% use a cell array for tps if you want to calculate std for different
% segments of the time course separately.
%


fnames = str2cell(fnames);

if ~iscell(tps)
    tps = {tps};
end


for i=1:length(fnames)
[ts,info] = BrikLoad(fnames{i});

      
for j=1:length(tps)
    
    if ~isempty(tps{j})
     ts_tmp = ts(:,:,:,tps{j});
    else
     ts_tmp = ts;
    end



std_ts = std(ts_tmp,0,4);
if ~exist('std2_sum','var')
    std2_sum = std_ts.^2;
    n = (std_ts>0);
else
    std2_sum = std_ts.^2+std2_sum;
    n = n+(std_ts>0);
end

end

end

brikdata = sqrt(std2_sum./n);
brikdata(isnan(brikdata)) = 0;

 nSb = 1;
info.DATASET_RANK(2) = nSb;  % the number of subbricks to save.
       info.BRICK_TYPES= 3*ones(1,nSb);  %1: short, 3: float
      
       info.BRICK_LABS = [];
       info.IDCODE_STRING = 'ave stdev';
       
       info.TYPESTRING = '3DIM_HEAD_ANAT';
       info.SCENE_DATA = [0 2 0];  % first 0: orig view; second 2: ANAT_EPI_TYPE; 0: matches TYPESTRING 
       
       info.BRICK_STATS=[];  %minimum and maximum values of the subbrick;
       
       info.BRICK_FLOAT_FACS = [];
       
       info.IDCODE_DATE = date;
       info.TAXIS_NUMS=[];
       history1 = 'tsstd(fnames,tps)';
       history2 = ['fnames = ',sprintf('%s\\n',fnames{:})];
       history3 = 'tps = ';
       for i=1:length(tps)
           history3 = [history3,num2str(tps{i})];
       end
       
       info.HISTORY_NOTE = [info.HISTORY_NOTE,'\n',history1,'\n',history2,'\n',history3];
       opt.OverWrite = 'y';

       opt.Prefix = out_name;   
       WriteBrik(brikdata,info,opt);







   

    