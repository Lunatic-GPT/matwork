function svm_block_ave(flist,t_trial,t_area_on,t_area_off,prefix)
% svm_block_ave(flist,t_trial,t_area_on,t_area_off,prefix)
% This function calculates the peak height and baseline values for each trial.
% flist: a cell array of afni time course file names.
% t_trial: number of data points for each trial.
% t_area_on: indices of data points used for peak height calculation.  The time indices are 1 based. 1 means the first time point after stimulus
%            onset.  use 0 for t_area_off is there is no off period (i.e. alternation
%            between two stimuli with no baseline periods) 
% t_area_off: indices of data points used for baseline calculation
% prefix: prefix of the output file.
% The number of bricks is 3 times the total number of blocks (n)
% the first n briks are areas of the peak.
% the second n briks are areas of the baseline.
% the third n briks are differences between the areas of the first and the
% second n blocks.

flist = str2cell(flist);

for i=1:length(flist)
    if i==1
        [data,info] = BrikLoad(flist{i});
    else
        tmp = BrikLoad(flist{i});
        data = cat(4,data,tmp);
    end
end

ntrials = size(data,4)/t_trial;
sz = size(data);

brik = zeros([sz(1:3),ntrials*3]); % peak areas, baseline areas, peak area -  baseline area

for i=1:ntrials
   
    brik(:,:,:,i) = mean(data(:,:,:,t_trial*(i-1)+t_area_on),4);
    if ~any(t_area_off==0) 
      brik(:,:,:,ntrials+i) = mean(data(:,:,:,t_trial*(i-1)+t_area_off),4);
    end
    brik(:,:,:,ntrials*2+i) = brik(:,:,:,i) - brik(:,:,:,ntrials+i); 
end

str = '';
for i=1:length(flist)
    str = [str,flist{i},' '];
end
history = sprintf('%s(%s,%s,%s,%s,%s)',mfilename,str,num2str(t_trial),num2str(t_area_on),num2str(t_area_off),prefix);

WriteBrikEZ(brik,info,history,prefix);


    