function d=Physio_Read(fname)

fid=fopen(fname,'r');
d=textscan(fid,'%d');

fclose(fid);


