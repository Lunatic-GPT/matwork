function convolve(varargin)
% convolve(ts_data1,stim_profile1[,ts_data2,stim_profile2 ...],prefix)
% the time series will corrected for offset before convolution.  The offset
% is calculated using the last four data points.

% 6/25/2009: remove ts_shift process.
      tic;
      for i=1:2:(nargin-2)
         [temp_sh,info] = BrikLoad(varargin{i});
         if i==1
             ts = zeros(size(temp_sh));
             len = size(ts,4);
         end
         
         mat_conv = reshape(varargin{i+1},[1,1,1,length(varargin{i+1})]);
         
         ts_conv = convn(temp_sh,mat_conv);
         ts = ts+ts_conv(:,:,:,1:len);
         
      end
      
      flist = sprintf('%s',varargin{1:2:(nargin-2)});
      stimulus = [];
      for i=2:2:nargin-1
          temp = ['stimulus ',num2str(i/2),': ', sprintf('%d ',varargin{i})];
          stimulus = [stimulus,temp,'\n'];
      end
      proc = ['processed in ', mfilename];
         info.BRICK_LABS = [];  
         info.HISTORY_NOTE = [info.HISTORY_NOTE, '\n',proc, '\ninput datasets:',flist,'\n',stimulus];
         opt.Prefix = varargin{end};
         opt.OverWrite = 'y';
         WriteBrik(ts,info,opt);
         disp([mfilename,' finish in ',num2str(toc),' s']);     
       
        
        