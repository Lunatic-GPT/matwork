function recon_MDEFT(d)

b=readbPar(fullfile(d,'method'),'PVM_EncSteps1',true);

a=readb_fid(d);

ra=a;
ra(:,b+length(b)/2+1,:)=a;

fra=fft2c(ra);

fra2=fft1c(fra,3);

order=readbPar(fullfile(d,'method'),'PVM_ObjOrderList');
fra3=fra2;
fra3(:,:,order+1,:)=fra2;

write_afni(abs(fra3),[d,'_recon']);
write_afni(abs(ra),[d,'_kspace']);
