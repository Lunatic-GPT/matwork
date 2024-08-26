function z = read_fid(fid_dir,format)
% z = read_fid(fname[,format])
% the output is a complex matrix of dimension np/2,ntraces,nblocks.
% where np is the number of data points for each fid.
if exist([fid_dir,'.fid'],'dir') 
    fid_dir = [fid_dir,'.fid'];
 end

 if ~isdir(fid_dir) && exist(fid_dir,'file')
    fname=fid_dir;
 else
    
 if exist(fullfile(fid_dir,'fid.orig'),'file')
  fname=fullfile(fid_dir,'fid.orig');
 else
   fname=fullfile(fid_dir,'fid');
 end
 end


fid = fopen(fname,'r','ieee-be');

fh = read_rawfheader(fid);

z = zeros(fh.np/2,fh.ntraces,fh.nblocks);

status = dec2bin(fh.status);

for i=1:fh.nblocks
  
    for j=1:fh.nbheaders
        
     bh(i) = read_rawbheader(fid);
    end
    
    if status(end-3) == '1' && fh.ebytes == 4 && ~exist('format','var')
      format = 'float';
    elseif status(end-3) == '0' && fh.ebytes == 2 && ~exist('format','var') 
      format = 'int16';
    elseif status(end-3) == '0' && fh.ebytes == 4 && ~exist('format','var')
      format = 'int32';
    end
    
    tmp = fread(fid,fh.np*fh.ntraces,format);
    z(:,:,i) = reshape(tmp(1:2:end-1),fh.np/2,fh.ntraces)+1i*reshape(tmp(2:2:end),fh.np/2,fh.ntraces); 
    
   % dr(:,:,i) = reshape(tmp(1:end/2),fh.np/2,fh.ntraces); 
    %di(:,:,i) = reshape(tmp(end/2+1:end),fh.np/2,fh.ntraces); 
    
end
fclose(fid);
