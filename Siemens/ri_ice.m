function img=ri_ice(dname,msize,nfiles)


if ~exist('msize','var') || isempty(msize)
    d=dir(fullfile(dname,'*.IceHead'));
    col=readxPar(fullfile(dname,d(1).name),'NoOfCols');
    msize(1)=str2num(col{1});
    
    row=readxPar(fullfile(dname,d(1).name),'NoOfRows');
    msize(2)=str2num(row{1});

end

dir_str=dir(fullfile(dname,'*.ima'));
if ~exist('nfiles','var')
    nfiles=length(dir_str);
end

for i=1:nfiles
    
 a=fopen(fullfile(dname,dir_str(i).name),'r');

d=fread(a,inf,'uint16');
d=reshape(d,msize);
img(:,:,i)=permute(d,[2,1,3,4]); % matrix column corresponds to COL in ICE if the file was saved with WriteToFile function

fclose(a);
end

if nargout==0
    save(dname,'img');
end
