function out = readimg(filename)

ind = findstr(filename,'.img');
if ~isempty(ind)
  filename = filename(1:ind-1);
end

ind = findstr(filename,'.hdr');
if ~isempty(ind)
  filename = filename(1:ind-1);
end

fmt=readhdr([filename,'.hdr'],'number format');

nd=readhdr([filename,'.hdr'],'number of dimensions');

mtx(1)=readhdr([filename,'.hdr'],'matrix size [1]');
mtx(2)=readhdr([filename,'.hdr'],'matrix size [2]');
mtx(3)=readhdr([filename,'.hdr'],'matrix size [3]');


% Open Data File
fid = fopen([filename,'.img'],'r','ieee-le');

if fid < 0
  error('Error opening data file');
end

switch fmt
  case 'short float'
    dtype = 'float32';
    otherwise
    error('Unsupported datatype');
end

% Read binary data
out = fread(fid,inf,dtype);

% Close data file
fclose(fid);

out=reshape(out,[mtx,length(out(:))/prod(mtx)]);

  
  
  