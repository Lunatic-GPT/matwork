function bheader = read_rawbheader(fid)

f = fread(fid,4,'int16');
bheader.scale = f(1);
bheader.status = f(2);
bheader.index = f(3);
bheader.mode = f(4);
bheader.ctcount = fread(fid,1,'int32');
a = fread(fid,4,'float');
bheader.lpval = a(1);
bheader.rpval = a(2);
bheader.lvl = a(3);
bheader.tlt = a(4);
