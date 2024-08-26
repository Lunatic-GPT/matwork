function d=read_phasefile(fname,np)
%d=read_fdf(fname[,f_afni])

fid = fopen(fname,'r','ieee-be');
dir_str=dir(fname);
d=fread(fid,dir_str.bytes/2,'float16');


nt=round(length(d)/np);

skip = length(d)-np*nt;
%skip=0;
r = d(1:end-skip);

d=reshape(r,[np,nt]);


