%function ase_fit(dname)

dname='EP2D_ASE_0011';



% images without z-shim
d=ri(dname,1);
extp(dname);
ns=readsPar([dname,'.pro'],'sSliceArray.lSize');


nro=size(d,1);
npe=size(d,2);
shimStep=readsPar([dname,'.pro'],'alFree[11]');
os=readsPar([dname,'.pro'],'adFree[12]');
nTE = size(d,4)/ns/os/shimStep;

d=reshape(d,[nro,npe,ns,shimStep*os,nTE]);

d2=squeeze(mean(d,4));

d2=interleave2linear(d2,3);



%% upload the first 36 images to server; run to3d; download the afni files;
%% number of TE

nte=readsPar([dname,'.pro'],'alFree[10]');
te=zeros(1,nte);

try
    te(1)=readsPar([dname,'.pro'],'alFree[13]');
catch
    te(1)=0;
end

for i=1:nte-1
    
    te(i+1)=readsPar([dname,'.pro'],sprintf('alFree[%d]',i+13));
    
end


%%
prefix='18_vr';

d=ri([prefix,'+orig']);

spatial_smoothing(d,0.5*[1,1,1],[prefix,'_sig0p5']);

d2=ri([prefix,'_sig0p5.mat']);

d2=reshape(d2,[64,64,36,3,6]);

d3=mean(d2,5);

[t2,M0]=T2map_LinearFit_data(d3(:,:,:,2:3),te(2:3));

cbv=1-d3(:,:,:,1)./M0;

r2=1./t2*1000;


save_nii(make_nii(r2),['R2pmap_',prefix,'_sig0p5.nii.gz']);

save_nii(make_nii(cbv),['CBV_',prefix,'_sig0p5.nii.gz']);
%% r2 map

prefix='17_zshim0_vr';
d=ri([prefix,'+orig']);

%spatial_smoothing(d,[0.5,0.5,0.5],[prefix,'_sig0p5']);

%d2=ri([prefix,'_sig0p5.mat']);
d2=d;
d2=reshape(d2,[64,64,36,3,6]);
d3=mean(d2,5);

[t2,M0]=T2map_LinearFit_data(d3(:,:,:,2:3),te(2:3));


cbv=1-d3(:,:,:,1)./M0;

r2=1./t2*1000;

%r2=interleave2linear(r2);
%cbv=interleave2linear(cbv);

save_nii(make_nii(r2),['R2pmap_',prefix,'_nosmooth.nii.gz']);

save_nii(make_nii(cbv),['CBV_',prefix,'_nosmooth.nii.gz']);