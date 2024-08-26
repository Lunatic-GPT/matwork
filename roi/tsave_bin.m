function tsave_bin(fnames,sti_file,trialTR,stimTR,out_name)

% tsave_bin(fnames,sti_file,trialTR,stimTR,out_name)
%
% fnames:  name of the input afni datasets
% sti_file: a file containing a matrix specifying the onset time of different trials
% trialTR: duration of one trial in units of TR.
% stimTR: number of time points during the on-period.
% out_name: prefix of the ave afni dateset.  The dataset names are
% [out_name,'_coh#+orig']

tic;
stim=load(sti_file);

if ~iscell(fnames)
    temp = fnames;
    fnames = cell(1);
    fnames{1} = temp;
end

ts = cell(1,length(fnames));

for i = 1:length(fnames)
    
    [ts{i},info] = BrikLoad(fnames{i});
 
end

nscan = length(fnames);
nTR = size(ts{1},4);
ncoh = size(stim,2);
nTrialCohExp = nscan*nTR/trialTR/ncoh; % number of trials of each coherence of the whole experiment


sz = size(ts{1});
bgTR = trialTR - stimTR;
nSb = 2*trialTR-stimTR;
  
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
       history1 = sprintf('tsave_bin(fnames,sti_file = %s,trialTR = %d,stimTR=%d,%s)',sti_file,trialTR,stimTR,out_name);
       history2 = ['fnames = ',sprintf('%s\\n',fnames{:})];
       info.HISTORY_NOTE = [info.HISTORY_NOTE,'\n',history1,'\n',history2];
       opt.OverWrite = 'y';
       


for icoh = 1:ncoh  % analyze separately for each coherence     
    noff1 = 0; % the number of 'off' periods after the coherent dots
    noff2= 0; % the number of 'off' periods before the coherent dots
    
    %find the time points when the current stimulus was on  
      ton = find(stim(:,icoh));  % time when the stimulus was on
      ntrial = length(ton)/stimTR;    % trials with this coherence
      ts_tmp = zeros([sz,ntrial]);
      
      for itrl = 1:ntrial
          ttemp = ton((itrl-1)*stimTR+1);
          iscan = ceil(ttemp/nTR);
          t_in = mod(ttemp,nTR);   % time within each scan.
              
          if(t_in+trialTR-1 <= nTR) && (t_in-bgTR >= 1)  % if there is a random dots period before and after the current coherent dots     
           noff1 = noff1+1;
           noff2 = noff2 +1;
           fst = t_in - bgTR;
           lst = t_in+trialTR-1;
           ts_tmp(:,:,:,:,itrl) = ts{iscan}(:,:,:,fst:lst);
          elseif t_in+trialTR-1 <= nTR %only after not before.
           fst = t_in;
           lst = fst + trialTR -1;
           noff1 = noff1+1;
           ts_tmp(:,:,:,bgTR+1:bgTR+trialTR,itrl) = ts{iscan}(:,:,:,fst:lst);
          elseif t_in-bgTR >= 1 % only before not after
              fst = t_in-bgTR;
              lst = t_in+stimTR-1;
              noff2 = noff2+1;
              ts_tmp(:,:,:,1:trialTR,itrl) = ts{iscan}(:,:,:,fst:lst);
          else
              error('This trial has no baseline condition'); 
          end         
      end
      ts_sum = sum(ts_tmp,5);
      brikdata(:,:,:,1:bgTR) = ts_sum(:,:,:,1:bgTR)/noff2;
      brikdata(:,:,:,bgTR+1:trialTR) = ts_sum(:,:,:,bgTR+1:trialTR)/ntrial;
      brikdata(:,:,:,trialTR+1:trialTR+bgTR) = ts_sum(:,:,:,trialTR+1:trialTR+bgTR)/noff1;
      
      opt.Prefix = [out_name,'_coh',num2str(icoh)];   
      WriteBrik(brikdata,info,opt);
         
end % loop through all the scans


%% save the average time series for each trial type.    

    disp([mfilename ' finish in ', num2str(toc), ' s']);