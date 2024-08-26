function [data,orig_delta] = reorient_data(data,orig_delta,orient)
%[data,orig_delta] = reorient_data(data,orig_delta,orient)
% permute data dimensions such that lr,ap and is are first, second and third
% dimensions in data, respectively.
if strcmp(orient,'sag') 
    ind = [3,2,1];
elseif strcmp(orient,'sag90')
    ind = [2,3,1];
elseif strcmp(orient,'cor')
    ind = [1,3,2];
elseif strcmp(orient,'cor90')
    ind = [3,1,2];
elseif strcmp(orient,'trans')
    ind = [1,2,3];
elseif strcmp(orient,'trans90') 
    ind = [1,2,3];
   % data = flipdim(data,2);
    %orig_delta(1,2) = orig_delta(1,2)+(size(data,2)-1)*orig_delta(2,2);
    %orig_delta(2,2)=-orig_delta(2,2);
else 
    ind = [1,2,3];
    %error('Unknown orientation');
end

ind2 = ind;
if ndims(data)==4
    ind2 = [ind,4];
end
    data = permute(data,ind2);
    if ~isempty(orig_delta)
    orig_delta = orig_delta(:,ind);
    end
    
    if strcmp(orient,'trans')
     
        data=flipdim(data,1);
    end