function res=myTextScan(fname)

fid=fopen(fname,'r');

a=textscan(fid,'%s');

res=a{1};

fclose(fid);