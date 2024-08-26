function recon_FieldMap(scand)

b=readbPar(fullfile(scand,'acqp'),'ACQ_dim');
if b~=3
    error('Not a 3D scan');
end


a=readb_fid(scand);

te=readbPar(fullfile(scand,'acqp'),'ACQ_echo_time');

a=reshape(a,[size(a,1),2,size(a,2),size(a,3)]);
a=permute(a,[1,3,4,2]);


of2=readbPar(fullfile(scand,'method'),'PVM_SPackArrPhase2Offset');
of1=readbPar(fullfile(scand,'method'),'PVM_SPackArrPhase1Offset');
of0=readbPar(fullfile(scand,'method'),'PVM_SPackArrReadOffset');
fov=readbPar(fullfile(scand,'method'),'PVM_Fov');


a=fov_shift_kspace(a,[-of0,-of1,-of2],fov);


fa=fft2c(a);
fa2=fft1c(fa,3);
fa2=flipdim(fa2,1);
fa2=flipdim(fa2,2);

write_afni(abs(fa2),[scand,'_recon']);
b=angle(conj(fa2(:,:,:,1)).*fa2(:,:,:,2))/diff(te)/2/pi*1000;  %b in hertz.
write_afni(b,[scand,'_B0']);

