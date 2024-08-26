function d=read_fdf_noheader(fname,dim,f_afni)


fid = fopen(fname,'r','ieee-be');
d=fread(fid,dim(1)*dim(2)*dim(3)*dim(4),'float32');
d = reshape(d,dim);

d=flipdim(d,2);
d=flipdim(d,3);
d = permute(d,[3,2,1,4]);

if exist('f_afni','var')
 [tmp,info]=BrikLoad(f_afni);
 prefix = fname(1:end-2);
 WriteBrikEZ(d,info,'read_fdf_noheader',prefix);
end

