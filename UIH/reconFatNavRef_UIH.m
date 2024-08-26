function reconFatNavRef_UIH(dname,nii_parent)
%% Three reference images will be constructed
% 1: grappa without replacement
% 2. grappa with replacement
% 3: use fully sampled data
% 
load([dname,'_smask.mat']);
load([dname,'_ref.mat']);

[~,suf]=strtok(dname,'_');
xmlname=name4pat(fullfile(dname,['*',suf,'.prot']),1);
xml=parseXML(xmlname);
lRefLine=readPar_uih(xml,'FatNavPATRefLines');
af=readPar_uih(xml,'FatNavPATFactor');

dacs=dref(:,end/2-lRefLine/2+1:end/2+lRefLine/2,end/2-lRefLine/2+1:end/2+lRefLine/2,:,1);
sz=size(dref);
im=zeros(sz(1),sz(2),sz(3),3,'single');



for j=1:size(dacs,1)
        
       
       tmp=shiftdim(dref(j,:,:,:),1);

       tmp=tmp.*squeeze(smask(1,:,:));

        kres=0*tmp;
       
        kres(inc_row,inc_col,:,:)=GRAPPA_new(tmp(inc_row,inc_col,:),squeeze(dacs(j,:,:,:)),[7,7],0.01,[af,af]);
     
        im(j,:,:,1) = sos(ifft1c(ifft1c(kres,1),2),3); %need to check before
        
        kres(end/2-lRefLine/2+1:end/2+lRefLine/2,end/2-lRefLine/2+1:end/2+lRefLine/2,:)=squeeze(dacs(j,:,:,:));

        im(j,:,:,2) = sos(ifft1c(ifft1c(kres,1),2),3); %need to check before


end

im(:,:,:,3)=sos(ifft1c(ifft1c(dref,2),3),4);


nii=load_untouch_niigz(nii_parent);
nii.img=im;
save_untouch_niigz(nii,[dname,'_ref.nii.gz']);







