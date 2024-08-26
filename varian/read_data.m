function d=read_data(fname,np)
%d=read_fdf(fname[,f_afni])

fid = fopen(fname,'r','ieee-be');
dir_str=dir(fname);
d=fread(fid,dir_str.bytes/4,'float32');


nt=round(length(d)/np/2);

skip = length(d)-np*nt*2;
%skip=0;
r = d(1:2:end-skip);
i=d(2:2:end-skip);
d=r+1i*i;
d=reshape(d,[np,nt]);


