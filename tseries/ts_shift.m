function data = ts_shift(data,npts)
%data = ts_shift(data,npts)
        % method 'end': shift the data so that the mean of the last npts data
        % points is equal to zero.
% the time is assumed to be the last dimension in data.
        sz = size(data);
        nd = length(sz);
        if nd ==4     
          offset = mean(data(:,:,:,end-npts+1:end),4);
          data = data - repmat(offset,[1,1,1,sz(4)]);
        elseif nd == 3
          offset = mean(data(:,:,end-npts+1:end),3);
          data = data - repmat(offset,[1,1,sz(3)]);
        elseif nd == 2
          offset = mean(data(:,end-npts+1:end),2);
          data = data - repmat(offset,[1,sz(2)]);
        else 
        disp([mfilename, 'time course not shifted for data dimensions other than 2,3 or 4']);
        end


