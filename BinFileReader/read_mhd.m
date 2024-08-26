function roi=read_mhd(fname)


fname=strtok2(fname,'.');

fid=fopen([fname,'.mhd'],'r');
list=textscan(fid,'%s');
fclose(fid);
pos=strmatch('DimSize',list{1});

if isempty(pos)
    error('Wrong mhd file format');
end

dim(1)=str2double(list{1}{pos+2});
dim(2)=str2double(list{1}{pos+3});
dim(3)=str2double(list{1}{pos+4});



f=fopen([fname,'.raw'],'rb');

d=fread(f,'*int16');

fclose(f);


  roi=reshape(d,dim);

roi=permute(roi,[2,1,3]);

%fname=strtok(fname,'.');
%save(fname, 'roi');


