function tsAreaHeight(flist,tp_area,scale,prefix)
%tsAreaHeight(flist,tp_area,scale,prefix)
% flist: list of files to compute the area.
% tp_area: time points to calculate area, 1 based.
% 6-3-2009 wrote it.
tic;
nts = length(flist);

data = cell(1,nts);
for i=1:nts
    [temp,info] = BrikLoad(flist{i});
    data{i} =  ts_shift(temp,'end');
end

sz = size(data{1});
area = zeros([sz(1:3),nts]);
height = zeros([sz(1:3),nts]);

for i=1:nts
    
    area(:,:,:,i) = mean(data{i}(:,:,:,tp_area),4)*scale;
    data_s = sort(data{i},4,'descend');
    height(:,:,:,i) = mean(data_s(:,:,:,1:3),4)*scale;
end

       info.DATASET_RANK(2) = nts;  % the number of subbricks to save.
       info.BRICK_TYPES=3*ones(1,nts);  %1: short, 3: float
       info.BRICK_LABS =[];
       info.BRICK_STATS =[];
       info.BRICK_FLOAT_FACS = zeros(1,nts);
       history = sprintf('tsAreaHeight(flist,tp_area,scale=%3.2f,prefix=%s)',scale,prefix);
       history3 = ['tp_area = ',num2str(tp_area)];
       history2 = ['flist =\n' sprintf('%s\\n',flist{:})];
       info.HISTORY_NOTE = [info.HISTORY_NOTE,'\n', history,'\n',history3,'\n',history2];
       Opt.OverWrite = 'y';
       Opt.Prefix = ['Height_',prefix];
       WriteBrik(height,info,Opt);
       
       Opt.Prefix = ['Area_',prefix];
       WriteBrik(area,info,Opt);

  disp([mfilename, ' finished in ',num2str(toc),' sec']);     
       

        