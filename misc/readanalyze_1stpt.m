function readanalyze_1stpt(prefix)
[a,d]=readanalyze(prefix);

writeanalyze(a(:,:,:,1),[prefix,'_1stpt'],d);
