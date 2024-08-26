function rotate_cbf_data(fname)

prefix=strtok(fname,'.');

[a,d]=readanalyze(prefix);
a2=permute(a,[2,1,3,4]);
d([2,1])=d([1,2]);

a2(a2<0)=0;
writeanalyze(a2,[prefix,'_rot'],d);

