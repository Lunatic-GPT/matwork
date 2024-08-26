function endian_l2b(fin,fout)
% endian_l2b(fin,fout)

fid=fopen(fin,'r','ieee-le');

a=dir(fin);
s=a.bytes/8;
d=fread(fid,s,'double');

fid2=fopen(fout,'w','ieee-be');
fwrite(fid2,d,'double');
fclose(fid);
fclose(fid2);