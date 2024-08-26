function a=rf_read(fname,v5)
%rf_read(fname,v5)

l =dir(fname);
npts = l.bytes/2;
f = fopen(fname,'r','ieee-be');
a=fread(f,npts,'int16');
fclose(f);

if exist('v5','var') && v5
    a = a(33:end-4);
end
figure;plot(a);

