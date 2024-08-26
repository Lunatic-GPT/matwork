root_dir = 'c:/users/xiaopeng/Desktop/cs_simulation';
load(fullfile(root_dir,'synthesize_kdata_mask2_nih_nois0.010_bold0.015.mat'));
ref=load('refg_0_10_30_6cy_TR2.0');
load(fullfile(root_dir,'mask_uniform_nc6.mat'));

m=reshape(m,[1,64,9,240]);
m=repmat(m,[64,1,1,1]);

FT2 = p2DFTxp(ones(size(m)));
im3=FT2'*z3;
z4=remove_mean_fullk(z3,m);
mz3=mean(z3,4);
tmp=z3-repmat(mz3,[1,1,1,240]);
tmp3=reshape(tmp,[64*64*9,240]);
tmp4=sqrt(mean(abs(tmp3).^2));

z5=mean_sampled_kt(z3,m);
z5=repmat(z5,[1,1,1,size(m,4)]);


im_init=zeros(size(m));
im_true=FT2'*z_true;
sl=4;
roi_ind=find(mask(:,:,sl,1)>0);

for i=1:length(roi_ind)
    [i1(i),i2(i)]=ind2sub([64,64],roi_ind(i));
end
ind_sparse=zeros(1,length(roi_ind));
for i=1:length(roi_ind)
  ind_sparse(i)=sub2ind([64,64,1,240],i1(i),i2(i),1,1);
end

%{
xfm=KLT(sl4);
FT2=p2DFTxp(ones(64,64,1,240));
res3=xfm*(FT2'*z3(:,:,sl,:));
res_true=xfm*(FT2'*z_true(:,:,sl,:));
%}


    norm1=sqrt(mean(abs(z3(:)-z_true(:)).^2));
    norm2=mean(abs(res3(ind_sparse)-res_true(ind_sparse))./abs(res_true(ind_sparse)));
    
    fprintf('l2 norm of error wrt noiseless image = %f\n',norm1);
    fprintf('mean fractional error at sparse indices = %f\n',norm2);
    

    
%sl4=run_cs((z3(:,:,sl,:)-z5(:,:,sl,:)),z3(:,:,sl,:)-z5(:,:,sl,:),z_true(:,:,sl,:)-z5(:,:,sl,:),ind_sparse,m(:,:,sl,:),xfm,0,0.01,10,false);
%sl4=run_cs(z3(:,:,sl,:),z3(:,:,sl,:),z_true(:,:,sl,:),ind_sparse,m(:,:,sl,:),xfm,0,0.1,10,false);


%sl4=ktFOCUSS_new(z3(:,:,sl,:)-z5(:,:,sl,:),z_true(:,:,sl,:)-z5(:,:,sl,:),ind_sparse,3,m(:,:,sl,:),0.5,0,20,2,1);

sl=4;
sl4=ktFOCUSS_new(z3(:,:,sl,:),z5(:,:,sl,:),ind_sparse,3,m(:,:,sl,:),0.5,0,20,2,1);

%ts=mean_roi(sl4,mask(:,:,sl,1));

sl4b=reshape(sl4,[64,64,1,240]);

im_fullk=im3;

tmp=im_true(:,:,sl,1);

thr=max(abs(tmp(:)))*0.05;

%% true vs recon
ts=mean_roi(abs(abs(64*sl4b./im_true(:,:,sl,:))-1),abs(tmp)>thr);
figure;plot(abs(ts));

%% true vs simul
ts=mean_roi(abs(abs(im_fullk(:,:,sl,:)./im_true(:,:,sl,:))-1),abs(tmp)>thr);
figure;plot(abs(ts));

figure;imshow(abs(abs(im_fullk(:,:,sl,30))-abs(im_true(:,:,sl,30))),[]);
figure;imshow(abs(64*abs(sl4b(:,:,1,30))-abs(im_true(:,:,sl,30))),[]);

%% noise correlation

a=fft(im_fullk(:,:,sl,:),[],4);
ma=mean_roi(abs(a),abs(tmp)>thr&mask(:,:,4,1)==0);
figure;plot(ma(2:end));


dim_fullk=im_fullk-repmat(mean(im_fullk,4),[1,1,1,240]);

a=repmat(dim_fullk(:,:,sl,2),[1,1,1,240]).*conj(dim_fullk(:,:,sl,:));
ma=mean_roi((a),abs(tmp)>thr&mask(:,:,4,1)==0);
figure;plot(abs(ma(2:end)));

a=fft(sl4b,[],4);
ma=mean_roi(abs(a),abs(tmp)>thr&mask(:,:,4,1)==0);
figure;plot(ma(2:end));

resid=sl4b*64-im_true(:,:,sl,:);
figure;imshow(abs(resid(:,:,10)),[]);

a=fft(sl4b*64-im_true(:,:,sl,:),[],4);
ma=mean_roi((a),mask(:,:,4,1)>0);
figure;plot(abs(ma(2:end)));

dsl4b=sl4b-repmat(mean(sl4b,4),[1,1,1,240]);
a=repmat(dsl4b(:,:,1,2),[1,1,1,240]).*conj(dsl4b(:,:,1,:));
ma=mean_roi(abs(a),abs(tmp)>thr&mask(:,:,4,1)==0);
figure;plot(abs(ma(3:end)));


%% klt transform
cz=zeros(size(z3));
cz(:,32,:,:)=z3(:,32,:,:);
sl4klt=KTFOCUSS_KLT_new(z3(:,:,sl,:),sl4,z_true(:,:,sl,:),ind_sparse,m(:,:,sl,:),0.5,0,20,2,1);

tsklt=mean_roi(sl4klt,mask(:,:,sl,1));
figure;plot(abs(tsklt));

%sl4=run_cs(z4(:,:,sl,:),im_init(:,:,sl,:),m(:,:,sl,:),xfm,0,0.01,false);

%tmp2=mean(z,4);



cc=zeros(64,64);
for i=1:64
    for j=1:64
     tmp=corrcoef(ref,abs(squeeze(sl4(i,j,:))));
     cc(i,j)=tmp(1,2);
    end
end


save('sparsity_random_nc_nois0.010_bold0.015','sl4','cc');
figure;imshow(cc(:,:,1),[]);