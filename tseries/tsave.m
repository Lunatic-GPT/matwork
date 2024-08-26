function tsave(fnames,sti_file,duration,npt_pre,out_name)

% tsave(fnames,sti_file,duration,npt_pre,out_name)
%
% fnames:  name of the input afni datasets
% sti_file: a file containing a matrix specifying the onset time of different trials
% duration: duration of one trial in units of TR.
% npt_pre: number of time points to include before each trial.
% out_name: prefix of the ave afni dateset.  The dataset names are
% [out_name,'_pat#+orig']
% 6/25/09:  remove detrend process. 

tic;
sti_seq = textread(sti_file,'','commentstyle','shell');
if sti_seq(1,end) ==0
    sti_seq(:,end) = [];
end
sti_seq = sti_seq';
sti_seq=sti_seq(:);


if ~iscell(fnames)
    temp = fnames;
    fnames = cell(1);
    fnames{1} = temp;
end

ts = cell(1,length(fnames));

for i = 1:length(fnames)
    
    [ts{i},info] = BrikLoad(fnames{i});
 
end

npt = size(ts{1},4);
npat = max(sti_seq);

sz = size(ts{1});
nSb = duration+npt_pre;
  
sz(4) = nSb;
  
      brikdata = zeros(sz);  

       info.DATASET_RANK(2) = nSb;  % the number of subbricks to save.
       info.BRICK_TYPES= 3*ones(1,nSb);  %1: short, 3: float
      
       info.BRICK_LABS = [];
       info.IDCODE_STRING = 'ave time series';
       
       info.TYPESTRING = '3DIM_HEAD_ANAT';
       info.SCENE_DATA = [0 2 0];  % first 0: orig view; second 2: ANAT_EPI_TYPE; 0: matches TYPESTRING 
       
       info.BRICK_STATS=[];  %minimum and maximum values of the subbrick;
       
       info.BRICK_FLOAT_FACS = [];
       
       info.IDCODE_DATE = date;
       info.TAXIS_NUMS=[];
       history1 = sprintf('tsave(fnames,%s,dur = %d,npt_pre=%d,%s)',sti_file,duration,npt_pre,out_name);
       history2 = ['fnames = ',sprintf('%s\\n',fnames{:})];
       info.HISTORY_NOTE = [info.HISTORY_NOTE,'\n',history1,'\n',history2];
       opt.OverWrite = 'y';
       
  for itr = 1:npat  %loop over the different trial types
      npre = 0;
      ind = find(sti_seq==itr);  
       
      nev = length(ind);
      ts_tmp = zeros([sz,nev]);
   
         
      for iev = 1:nev  %loop over all the events of the same trial
                ifl = floor((ind(iev)*duration-1)/npt)+1; % determine which file   
                fst = (ind(iev)-1)*duration+1-(ifl-1)*npt;
                lst = fst+duration-1; 
                     
                if fst >npt_pre % not the first trial in the scan.
                     fst = fst-npt_pre;
                     npre=npre+1;
                     ts_tmp(:,:,:,:,iev) = ts{ifl}(:,:,:,fst:lst);
                else
                    ts_tmp(:,:,:,npt_pre+1:npt_pre+duration,iev) = ts{ifl}(:,:,:,fst:lst);
                end   
      end
           
         ts_sum = sum(ts_tmp,5);
         brikdata(:,:,:,1:npt_pre) = ts_sum(:,:,:,1:npt_pre)/npre;
         brikdata(:,:,:,1+npt_pre:duration+npt_pre) = ts_sum(:,:,:,1+npt_pre:duration+npt_pre)/nev;
                      
         opt.Prefix = [out_name,'_pat',num2str(itr)];   
         WriteBrik(brikdata,info,opt);
         
  end  
       
     
%% save the average time series for each trial type.    

    disp([mfilename ' finish in ', num2str(toc), ' s']);