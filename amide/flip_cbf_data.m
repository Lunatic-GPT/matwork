function flip_cbf_data(fname)

prefix=strtok(fname,'.');

a=readanalyze(prefix);
a2=flipdim(a,2);
writeanalyze(a2,[prefix,'_flip'],[0.3333,0.3333,1]);

