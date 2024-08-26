function fheader = read_rawfheader(fid)

f = fread(fid,6,'int32');
fheader.nblocks = f(1);
fheader.ntraces = f(2);
fheader.np = f(3);
fheader.ebytes = f(4);
fheader.tbytes = f(5);
fheader.bbytes = f(6);

a=fread(fid,2,'int16');

fheader.vers_id = a(1);
fheader.status = a(2);
fheader.nbheaders = fread(fid,1,'int32');

