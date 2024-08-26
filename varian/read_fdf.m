function d=read_fdf(fname,f_afni)
%d=read_fdf(fname[,f_afni])

fid = fopen(fname,'r','ieee-le');
    a=textscan(fid,'%s%s%s%s%s%s',10);
    i1=a{4}{7}(2:end-1);
   
    i1=str2double(i1);
    if isempty(a{6}{7})
         i2=a{5}{7}(1:end-2);
         
    i2=str2double(i2);
        i3=1;
    else
         i2=a{5}{7}(1:end-1);
      i3=a{6}{7}(1:end-2);
      
    i2=str2double(i2);
    i3=str2double(i3);
    end
    


fseek(fid,-i1*i2*i3*4,1);
d=fread(fid,i1*i2*i3,'float32');
d = reshape(d,[i1,i2,i3]);
d=flipdim(d,2);
d=flipdim(d,3);

d = permute(d,[3,2,1]);
fclose(fid);
if exist('f_afni','var')
 [tmp,info]=BrikLoad(f_afni);
 prefix = fname(1:end-4);
 WriteBrikEZ(d,info,'read_fdf',prefix);
end

