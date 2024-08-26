function tsAlign(p)
% tsAlign(p)
%
% fnames:  name of the input afni datasets
% callback function of tsAlign_gui
% 6/25/09:  remove detrend process. 
% 6/13/09: modified to allow variable durations of different trial types.
%          different trial types are now saved in different files.
% 10/06/2010: allow baseline to be calculated in each trial and remove 
% linear drifts in each trial.  To use these feature, detrended and normalized time
% courses should be used. 
% 04-16-11: allow unequal number of trials across scans.
tic;
%set(gcbo,'Enable','off');
fnames_str = get(p,'time series');
  temp = textscan(fnames_str,'%s','Delimiter',',');
  fnames = temp{1};

brik = get(p,'subbrik range');  
sti_file = get(p,'stimulus sequence file');
duration = get(p,'trial duration (TR)');
bl_dur = get(p,'baseline duration (TR)');  % 0 or empty to skip adjustments.

out_name = get(p,'prefix for output files');

sti_seq = textread(sti_file,'','commentstyle','shell');
if sti_seq(1,end) ==0
    sti_seq(:,end) = [];
end

fnames = str2cell(fnames);


if size(sti_seq,1) == 1 || size(sti_seq,2) == 1
     sti_seq = reshape(sti_seq,length(sti_seq)/length(fnames),length(fnames));
     sti_seq = sti_seq';
end

if length(fnames) ~=size(sti_seq,1)
    error('Number of data files does not match that in the sequence file');
end


dim3 = (~isempty(strfind(fnames{1},'+orig')));
  
ts = cell(1,length(fnames));

for i = 1:length(fnames)
    if dim3
        if ~isempty(brik)
          [ts{i},info] = BrikLoadf([fnames{i},'[',brik,']']);
        else
           [ts{i},info] = BrikLoadf(fnames{i});
        end
    else
      tmp = load(fnames{i});
      ts{i} = shiftdim(tmp',-2);
    end
end



ntrl = max(sti_seq(:));

if length(duration) ==1
    duration = duration*ones(1,ntrl);
end


for itrl =1:ntrl
  
  brikdata = [];  
  
    
  for i=1:size(sti_seq,1)
      t_ind = 0;
      
    for j=1:size(sti_seq,2)
                  
        if sti_seq(i,j)==0
            continue;
        end
        if itrl ~= sti_seq(i,j)
            t_ind=t_ind+duration(sti_seq(i,j));
            continue;
        end
        
        if t_ind == 0
            ts_tmp = baseline_shift(ts{i}(:,:,:,t_ind+1:t_ind+duration(itrl)),bl_dur);
        else
            ts_tmp = baseline_detrend(ts{i}(:,:,:,t_ind+1-bl_dur:t_ind+duration(itrl)),bl_dur);
        end
        
        if isempty(brikdata)    
            brikdata = ts_tmp;
        else
          brikdata(:,:,:,end+1:end+duration(itrl)) = ts_tmp;
        end
        t_ind = t_ind + duration(itrl);
        disp([i,j]);
    end
    if t_ind ~= length(ts{i})
        error('time series length error');
    end
    
 end

if dim3
       info.BRICK_LABS = [];
       info.IDCODE_STRING = 'aligned time series';
       
       info.TYPESTRING = '3DIM_HEAD_ANAT';
       info.SCENE_DATA = [0 2 0];  % first 0: orig view; second 2: ANAT_EPI_TYPE; 0: matches TYPESTRING 
       
       info.BRICK_STATS=[];  %minimum and maximum values of the subbrick;
       
       info.BRICK_FLOAT_FACS = [];
       
       info.IDCODE_DATE = date;
       info.TAXIS_NUMS=[];
       history1 = sprintf('tsAlign(fnames,%s,durations,%s)',sti_file,out_name);
       history2 = ['fnames = ',sprintf('%s\\n',fnames{:})];
       history3 = ['duration = ',sprintf('%d,',duration)];
       info.HISTORY_NOTE = [info.HISTORY_NOTE,'\n',history1,'\n',history2,'\n',history3];

       disp('saving the aligned time series');
       nSb = size(brikdata,4);
       info.DATASET_RANK(2) = nSb;  % the number of subbricks to save.
       info.BRICK_TYPES= 3*ones(1,nSb);  %1: short, 3: float
      
       
       opt.Prefix = [out_name,'_p',num2str(itrl)];
       opt.OverWrite = 'y';
       
       WriteBrik(brikdata,info,opt);
     
else
      tmp = squeeze(brikdata);
      if size(tmp,2) ~= 1
          tmp = tmp';
      end
      save_mat_int(tmp,sprintf('%s_p%d.1D',out_name,itrl));
end

end
 % set(gcbo,'Enable','on');
     disp([mfilename ' finish in ', num2str(toc), ' s']);
         
         
function ts_tmp = baseline_shift(ts,bl)
% the length of ts is equal to trial duration

if isempty(bl) || bl == 0;
    ts_tmp = ts;
    return;
end

len = size(ts,4);
ts = ts+1; %convert back to original signal intensity in units of mean scan intensity
bl_trl = mean(ts(:,:,:,end-bl+1:end),4); %calculate trial mean

ts_tmp = ts./repmat(bl_trl,[1,1,1,len])-1;
ts_tmp(isnan(ts_tmp)) = 0;


function ts_tmp = baseline_detrend(data,bl_dur)

% the length of ts is equal to trial duration plus bl.  
% the first bl points are from the previous trial

if isempty(bl_dur) || bl_dur == 0;
    ts_tmp = data;
    return;
end

sz = size(data);
len = size(data,4) - bl_dur;
tps = [1:bl_dur,len+1:len+bl_dur];
order = 1;
 for i=1:size(data,1)
            for j=1:size(data,2)
                for k=1:size(data,3)
                    
                    ts = squeeze(data(i,j,k,tps));
                    p = polyfit(tps,ts',order);
                    t = 1:sz(4);
                    
                    bl = zeros(1,sz(4));
                    for ip=0:order
                       bl = bl+t.^(order-ip)*p(ip+1);
                    end
                    bl = shiftdim(bl,-2);                       
                    data(i,j,k,t) = data(i,j,k,t) - bl + squeeze(mean(data(i,j,k,len+1:len+bl_dur),4)); 
                end
            end
         %   disp(i);
 end
 
 
 ts_tmp = baseline_shift(data(:,:,:,bl_dur+1:end),bl_dur);
        
 


            
            
        
       
  