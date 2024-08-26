roi=load('mask_highK_slice1.mat','d');
roi=roi.d(:,:,1);

tmp=reshape(dsb2,[458,540,32*3]);

sd=zeros(32,3);
for i=1:32*3
    tmp2=tmp(:,:,i);
    sd(i)=std(real(tmp2(roi>0)));
end

%%
tmp=sos(dsb2,3);
[tmp2,ind_max]=max(tmp(:));
ind_max=ind2subb(size(tmp),ind_max);

%ind_max=[141,273];
lc=80;  % calibration lines
iread=ind_max(1)-lc/2:ind_max(1)+lc/2-1;


Ks=-2:2;
Ksy=-2:2;
iph=ind_max(2)-lc/2:ind_max(2)+lc/2-1;
kernel_all=[];
for i=1:100
    data_test=dsb2(iread,iph,:,:,1)+randn_white([length(iread),length(iph),32,3]);
    
    kernel_tmp=sliceGRAPPAKernel(data_test,Ks,Ksy);
    
    kernel_tmp=reshape(kernel_tmp,[32,length(Ks),length(Ksy),32,3]);
    kernel_all(:,:,:,:,i)=squeeze(kernel(29,:,:,:,:));
    disp(i);
end

    
%%
[tmp,ind_max]=max(kernel(:));
ind=ind2subb(size(kernel),ind_max);
    