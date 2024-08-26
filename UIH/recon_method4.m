function img=recon_method4(fd,motionPar,shotIndex,FOV,ro_pe_par,ds)
%fd should be in the k-space for all dimensions

voxsize=FOV./size(fd(:,:,:,1));
mask=squeeze(sos(sos(fd,1),4))>1;
[rows,cols]=fullySampledRegion(mask);

mtmp = 0*mask;
mtmp(rows,cols)=1;
mtmp = logical(mtmp(:));

[k,data]=get_k_data_moco_shotIndex_zbb(fd,motionPar,shotIndex,FOV,ro_pe_par,ds);
sel=sos(sos(data,3),1)>0;

k2=squeeze(k(:,sel,:));

%k2=k(:,sel,abs(ro_pe_par));


data2=data(:,sel,:);

calib=squeeze(data(:,mtmp,:));
kc=squeeze(k(:,mtmp,:));

% [kc,calib]=get_k_data_moco_shotIndex_zbb(fd(:,rows,cols,:),motionPar,shotIndex(rows,cols),FOV,ro_pe_par,ds);
% sel=sos(sos(calib,3),1)>0;
% kc2=squeeze(kc(:,sel,:));
% calib2=calib(:,sel,:);
    
%%
sz=size(fd);
img=zeros(sz);

par.ksize=[6,6];
par.wnthresh = 1.8;
par.eigThresh_im = 0.9;
par.nIterCG=15;
par.nIterSAKE=100;
par.imSize = [sz(2),sz(3)];
par.calSize = [sum(rows),sum(cols)];


kc_c = squeeze(kc(end/2+1,:,2:3));
dc_c = squeeze(calib(end/2+1,:,:));
ESP=get_ESP_SAKE_nufft(dc_c,kc_c,par);

ltime=tic;
%parpool(8);
for i=1:sz(1)

    % k space coordinates
    ki = squeeze(k2(i,:,2:3));
    % k space data
    kdi = squeeze(data2(i,:,:));

    kci = squeeze(kc(i,:,2:3));
    kcdi = squeeze(calib(i,:,:));
    [~,~,img(i,:,:,:)] = do_SAKE_ESPIRiT_NUFFT(ki,kdi,kci,kcdi,par,ESP);

    time_left(i,sz(1),toc(ltime));
    %         [res,resw,img(i,:,:,:)] = do_SAKE_ESPIRiT_NUFFT(ki,kdata,kci,kcdata,par);
end
img = ifft1c(img,1);
save recon_method4_img img  -v7.3;
img = sos(img,4);
end