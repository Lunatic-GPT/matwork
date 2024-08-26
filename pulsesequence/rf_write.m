function rf_write(fname,data)
%rf_read(fname,data)

f = fopen(fname,'w','ieee-be');
fwrite(f,data,'int16');
fclose(f);
