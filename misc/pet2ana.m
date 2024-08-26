function pet2ana(fname,brik)
% brik 1 based.
prefix=strtok(fname,'.');
fname=[prefix,'.hdr'];

en=readhdr(fname,'imagedata byte order');

fmt=readhdr(fname,'number format');

nbytes=readhdr(fname,'number of bytes per pixel');

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

        
sz(1)=readhdr(fname,'matrix size [1]');
sz(2)=readhdr(fname,'matrix size [2]');
sz(3)=readhdr(fname,'matrix size [3]');

res(1)=readhdr(fname,'(mm/pixel) [1]');
res(2)=readhdr(fname,'(mm/pixel) [2]');
res(3)=readhdr(fname,'(mm/pixel) [3]');
res(4)=1000;
sz(4)=readhdr(fname,'number of time frames');


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
d=flipdim(d,2);
%d = permute(d,[2,1,3,4]);
if ~exist('brik','var')
writeanalyze(d,[prefix,'_ana'],res,fmt);
else
    if isa(brik,'char')
    brik=str2num(brik);
    end
    if any(brik>sz(4))
      brik=sz(4);
    end
writeanalyze(sum(d(:,:,:,brik),4),sprintf('%s_t%d_ana',prefix,brik),res,fmt);

fprintf('save #%d of %d volumes\n',brik,size(d,4));
end    
