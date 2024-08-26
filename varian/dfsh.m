function dfsh(fid_prefix,lsfid)
%dfsh(fid_prefix,lsfid)
a=read_fid([fid_prefix,'.fid/fid']);
a=abs(a);
a=squeeze(a(lsfid+1:end,:,:));
figure;plot(a(:));

