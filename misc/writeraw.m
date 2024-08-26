function writeraw(d,fname,precision)
fid=fopen(fname,'w');

if exist('precision','var')
 fwrite(fid,d(:),precision);
else
  fwrite(fid,d(:));
end
    
fclose(fid);

