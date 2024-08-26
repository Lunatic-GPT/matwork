function recon_FLASH3D(scand)

b=readbPar(fullfile(scand,'acqp'),'ACQ_dim');
if b~=3
    error('Not a 3D scan');
end


a=readb_fid(scand);

sz=size(a);
NI=readbPar(fullfile(scand,'acqp'),'NI');
a=reshape(a,[sz(1),NI,sz(2:3)]);

a=permute(a,[1,3,4,2]);
%a=reshape(a,[sz(1),sz(3),sz(2),sz(4)]);

%a=permute(a,[1,3,2,4]);

fa=fft2c(a);
fa2=fft1c(fa,3);

write_afni(abs(fa2),[scand,'_recon']);
