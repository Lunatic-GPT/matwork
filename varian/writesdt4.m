function writesdt4(a,filename)
% writesdt4(a,filename)
sdtfile=strcat(filename,'.sdt');
fid=fopen(sdtfile,'w','ieee-le');
fwrite(fid,a,'float32');
fclose(fid);

sprfile=strcat(filename,'.spr');
fid=fopen(sprfile,'w');
numDimStr=ndims(a);
fprintf(fid,'numDim:%s\n',num2str(numDimStr));
dimStr=size(a);
fprintf(fid,'dim:%s\n',num2str(dimStr));
fprintf(fid,'dataType: REAL\n');
fprintf(fid,'fidName: %s\n',strcat(filename,'.sdt'));
fprintf(fid,'sdtOrient: ax\n');
fprintf(fid,'endian: ieee-le\n');
%fprintf(fid,'Real2WordScale: 1\n');

fclose(fid);

