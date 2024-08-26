function [ts_mean,ts_std]=tsroieb(xdata,ts_fnames,mask,ts_scale,barplot,out_fname)
% tsroi2eb(xdata,ts_fnames,mask[,ts_scale,barplot,out_fname])
% ts_fnames: a cell array of dataset names.   
% ts_scale: an array to scale each dataset.
tic;

if ~exist('barplot','var')
    barplot = 0;
end
if ~iscell(ts_fnames)
    temp = ts_fnames;
    ts_fnames =cell(1);
    ts_fnames{1}= temp;
end

nf = length(ts_fnames);

if  ~exist('ts_scale','var')|| isempty(ts_scale) 
    ts_scale =ones(1,nf);
end

fsize = 14;
for fi =1:nf
 
    cmd = sprintf('3dmaskave -mask ''%s'' ''%s'' > maskave_temp.1D',mask,ts_fnames{fi});
    unix(cmd);
    [ts_tmp,nvox_temp,temp_str] = textread('maskave_temp.1D','%f [%d %s');
    if ~strcmp(temp_str,'voxels]')
        error('3dmaskave output file error');
    end
    ts(:,fi) = ts_tmp(:)*ts_scale(fi);
end


ts_len = length(xdata);
ntrials = size(ts,1)/ts_len;
ts = reshape(ts,[ts_len,ntrials,fi]);
ts_mean = squeeze(mean(ts,2));
ts_std = squeeze(std(ts,0,2))/sqrt(ntrials);

ts_mean = ts_shift(ts_mean);
nvox = nvox_temp(1);
 % figure;
  if barplot
      h = barweb(ts_mean,ts_std);
      set(h,'BarWidth',1);
      colormap([0.6,0.6,0.6;0,0,1;0.6,0,0]);
  else
       
      % symb = {'k+','bs','r*','g>'};
       % symb = {'ks','b*','r+','g>','m','c'};
        symb = {'k','b','r','g','m','c'};    
        %symb = {'k-','b*','rs','k-','k-'};
       for i=1:size(ts_mean,2)
        errorbar(xdata,ts_mean(:,i),ts_std(:,i),symb{i});
        hold on;
       end
  end
  
  legend('toggle');
  down = min(ts(:));
  up = max(ts(:));
  nvox_str =sprintf('n=%d',nvox);
 text(ts_len/2,down+(up-down)/4,nvox_str,'FontSize',fsize);
 mask = strrep(mask,'_','\_');
 title(mask);
 set(gca,'box','off');
 
 %xlim([0,15]);
 %ylim([0,up*1.1]);
if exist('out_fname','var')
  save(out_fname,'ts','mask','ts_fnames','nvox');
end

       
disp([mfilename ' finish in ', num2str(toc), ' s']);
     
    

function ts = ts_shift(ts)

   offset = mean(ts(end-2:end,:),1);
    
   ts = ts - repmat(offset,[size(ts,1),1]);
      
         
         
         
         


