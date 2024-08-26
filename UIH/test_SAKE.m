load  8_fse_mx_3d_fatnav2_ACS5_fd170.mat

[rows,cols]=fullySampledRegion(sos(fd_170,3)>0);

kSize=[6,6];
 tmp = im2row(fd_170(rows,cols,:),kSize); [tsx,tsy,tsz] = size(tmp);
Nc=size(fd_170,3);
wnthresh=2;

    A = reshape(tmp,tsx,tsy*tsz);
    
    % SVD thresholding
    [U,S,V] = svd(A,'econ');
   
  

figure;


subplot(211);plot(1:size(S,2),diag(S),'LineWidth',2);
hold on, 
plot([wnthresh,wnthresh]*prod(kSize),[0,S(1)],'g--','LineWidth',2)
legend('singular vector value','threshold')
title('');
subplot(212);
imagesc(abs(V));


