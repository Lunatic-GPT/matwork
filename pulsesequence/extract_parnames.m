function par = extract_parnames(fname)
% extract the parameter names enclosed in "" in a varian pulse sequence. 
fid = fopen(fname,'r');
par = {};
a=textscan(fid,'%s');
ar=a{1};
for i=1:length(ar)
  ind = strfind(ar{i},'"');
  
  if ~isempty(ind)
    if length(ind)~=2
        error('file format error');
    end
    
    par{end+1} = ar{i}(ind(1)+1:ind(2)-1);
  end
end
	   