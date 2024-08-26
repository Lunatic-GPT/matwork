function [position,normal,inplanerot]=get_slicePosition(fname)

fid=fopen(fname,'r');
position=fscanf(fid,'position=(%f,%f,%f)\n');
normal=fscanf(fid,'normal=(%f,%f,%f)\n');
inplanerot=fscanf(fid,'inplanerotation=(%f)\n');


fclose(fid);