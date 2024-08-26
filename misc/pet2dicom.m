function pet2dicom(slice_combine)

if ~exist('slice_combine','var')
    slice_combine=false;
end

dir_str=dir('*.hdr');


for i=1:length(dir_str)
    
    fname=dir_str(i).name;
    prefix=fname(1:end-4);
en=readhdr(fname,'imagedata byte order');

fmt=readhdr(fname,'number format');

nbytes=readhdr(fname,'number of bytes per pixel');
nbytes=str2num(nbytes);

switch fmt
    case 'short float'
     fmt = 'float';
    case 'signed integer'
     if nbytes==2
     fmt = 'short';
     else
     fmt='int';
     end
    otherwise
        error('unknown format');
end

         


sz(1)=str2num(readhdr(fname,'matrix size [1]'));
sz(2)=str2num(readhdr(fname,'matrix size [2]'));
sz(3)=str2num(readhdr(fname,'matrix size [3]'));

res(1)=str2num(readhdr(fname,'(mm/pixel) [1]'));
res(2)=str2num(readhdr(fname,'(mm/pixel) [2]'));
res(3)=str2num(readhdr(fname,'(mm/pixel) [3]'));
sz(4)=str2num(readhdr(fname,'number of time frames'));

info.StudyTime = readhdr(fname,'study time');
info.StudyDate=readhdr(fname,'study date');
info.StudyDate = strrep(info.StudyDate,'/','.');
    info.StudyID=readhdr(fname,'study id');

info.ImageDate=info.StudyDate;
info.PatientName=readhdr(fname,'patient name');
info.PatientID=readhdr(fname,'patient ID');
info.PatientWeight=readhdr(fname,'patient weight');
info.PatientOrientation=readhdr(fname,'patient orientation');
info.PatientOrientation=upper(strrep(info.PatientOrientation,'_','\'));

info.SliceOrientation=readhdr(fname,'slice orientation');

str=readhdr(fname,'scaling factor (mm/pixel) [1]');

info.NumberOfFrames=sz(4);
info.PixelSpacing=res(1:2);
if strcmp(en,'LITTLEENDIAN') ||  strcmp(en,'littleendian')
    fid = fopen([prefix,'.img'],'r','l');
elseif strcmp(en,'BIGENDIAN') || strcmp(en,'bigendian')
    fid = fopen([prefix,'.img'],'r','b');
else
    error('unknown endian');
end

d = fread(fid,fmt);

fclose(fid);

d=reshape(d,sz);
d = permute(d,[2,1,3,4]);


if ~slice_combine
    for k=1:size(d,3)
     status=    dicomwrite(d(:,:,k,:),sprintf('%s-%03d.dcm',prefix,k),info,'SOPClassUID','1.2.840.10008.5.1.4.1.1.128', 'CreateMode', 'copy');        
     disp(sprintf('%s: %d frame(s), slice %d.',prefix,size(d,4),k));
    end
else
ncol=sqrt(size(d,3));
ncol=ceil(ncol);
nrow=ceil(size(d,3)/ncol);

sz=size(d);

 for j=1:size(d,4)
    

    d2=zeros(sz(1)*nrow,sz(2)*ncol);
    for k=1:size(d,3)
       col=mod(k-1,ncol)+1;
       row=ceil(k/nrow);
     d2(sz(1)*(row-1)+1:sz(1)*row,sz(2)*(col-1)+1:sz(2)*col)=d(:,:,k,j);
    end
   dicomwrite(d2,sprintf('%s_%03d.dcm',prefix,j));  
   disp(sprintf('%s_%03d.dcm.',prefix,j));
 end

end

end
