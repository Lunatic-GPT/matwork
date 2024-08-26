function tsHeight(flist,pos,prefix)
%tsHeight(flist,pos,prefix)
% flist: list of files to compute the height.
%pos: 0: use maximum value otherwise the time point at pos and neighbouring
%two time points to calculate the height. 1 based.
% 6-15-2009 wrote it.
tic;
nts = length(flist);

data = cell(1,nts);
for i=1:nts
    [temp,info] = BrikLoad(flist{i});
    data{i} =  ts_shift(temp,'end');
end

sz = size(data{1});
height = zeros([sz(1:3),nts]);

for i=1:nts
    if pos == 0
      data_s = sort(data{i},4,'descend');
      height(:,:,:,i) = mean(data_s(:,:,:,1:3),4);
    else
      height(:,:,:,i) = mean(data_s(:,:,:,pos-1:pos+1),4);  
    end    
end

       info.DATASET_RANK(2) = nts;  % the number of subbricks to save.
       info.BRICK_TYPES=3*ones(1,nts);  %1: short, 3: float
       info.BRICK_LABS =[];
       info.BRICK_STATS =[];
       info.BRICK_FLOAT_FACS = zeros(1,nts);
       history = sprintf('tsHeight(flist,prefix=%s)',prefix);
       history2 = ['flist =\n' sprintf('%s\\n',flist{:})];
       info.HISTORY_NOTE = [info.HISTORY_NOTE,'\n', history,'\n',history,'\n',history2];
       Opt.OverWrite = 'y';
       Opt.Prefix = prefix;
       WriteBrik(height,info,Opt);
       
  disp([mfilename, ' finished in ',num2str(toc),' sec']);     
       
function data = ts_shift(data,method)
        % method 'start': shift the data so that the first data point is equal to
        % 0.
        % method 'end': shift the data so that the mean of the last 4 data
        % points is equal to zero.
        
        sz = size(data);
        if length(sz) ~=4
         error([mfilename, ' only works for 3 dimensional data']);
        end
     
        switch method
            case 'sta'
             offset = data(:,:,:,1);
            case 'end'
             offset = mean(data(:,:,:,end-3:end),4);
            otherwise
             error([mfilename,' unknown method']);   
        end

        data = data - repmat(offset,[1,1,1,sz(4)]);
        