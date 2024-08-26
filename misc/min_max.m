function a = min_max(data,dim)
%a = min_max(data) 
% a is a 1 by 2 row vector containing min and max values of data.


 if ~exist('dim','var')
  tmp = data(:);
  a = zeros(1,2);
  a(1) = min(tmp);
  a(2) = max(tmp);
 else
  amin = min(data,[],dim);
  amax = max(data,[],dim);
  a=cat(dim,amin,amax);
 end