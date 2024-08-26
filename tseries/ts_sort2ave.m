function ts_sort2ave(fordered,nTR_trial,add_points,prefix,exclude)
% ts_mn = ts_sort2ave(fordered,nTR_trial,add_points,prefix,exclude)
% fordered:  fname of the ordered time series. Can have more than one files
% exclude: the trials to be excluded, 1 based;
% prefix: set to empty if do not generate new files.
% add_points: add add_points time points at the beginning of each scan.
tic;
fordered=str2cell(fordered);

for i=1:length(fordered)
    
    if ~any(fordered{i}=='+')
        fordered{i}=[fordered{i},'+orig'];
    end
    
     [tmp,info] = BrikLoadf(fordered{i});
     if add_points>0
         tmp = cat(4,tmp(:,:,:,1:add_points),tmp);
     end
     if i==1
      data = tmp;
     else
      data = cat(4,data,tmp);
    end
end

nTR = size(data,4);


if ~exist('exclude','var')
    exclude=[];
end

if mod(nTR,nTR_trial) ~= 0 
    error('mod(nTR,nTR_trial) ~= 0 ');
end

ntrl = nTR/nTR_trial;
sz=size(data);
data2 = reshape(data,[sz(1:3),nTR_trial,ntrl]);

include = setdiff(1:ntrl,exclude);

brik_data = mean(data2(:,:,:,:,include),5);
           
files=sprintf('%s;',fordered{:});
history = sprintf('ts_sort2ave(%s,%d,addpoints=%d,%s,%s)',files,nTR_trial,add_points,prefix,num2str(exclude));
WriteBrikEZ(brik_data,info,history,prefix);

disp([mfilename ' finish in ', num2str(toc), ' s']);


