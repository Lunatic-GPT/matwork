function save_slicePosition(fpos,pos,norm,rotAngle)
%save_slicePosition(fpos,pos,mat)

fid=fopen(fpos,'w');
fprintf(fid,'position=(%f,%f,%f)\n',pos);
fprintf(fid,'normal=(%f,%f,%f)\n',norm);

fprintf(fid,'inplanerotation=(%f)\n',rotAngle); % angle in radian
fclose(fid);

