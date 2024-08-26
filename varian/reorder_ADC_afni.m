function reorder_ADC_afni(prefix)
[a,info]=BrikLoad([prefix,'+orig']);
b(:,:,:,1)=mean(a(:,:,:,1:7:56),4);
b(:,:,:,2)=mean(a(:,:,:,2:7:56),4);
b(:,:,:,3)=mean(a(:,:,:,3:7:56),4);
b(:,:,:,4)=mean(a(:,:,:,4:7:56),4);
b(:,:,:,5)=mean(a(:,:,:,5:7:56),4);
b(:,:,:,6)=mean(a(:,:,:,6:7:56),4);
b(:,:,:,7)=mean(a(:,:,:,7:7:49),4);
img(:,:,:,1)=b(:,:,:,7);
img(:,:,:,2)=mean(b(:,:,:,1:6),4);
img(:,:,:,3)=img(:,:,:,2)./img(:,:,:,1);

dim1=size(img,1);
dim2=size(img,2);
dim3=size(img,3);
c_tmp=img(:,:,:,1);
SI_thres=0.05*max(c_tmp(:));
b=[5,1200];
tmp=zeros(dim1,dim2,dim3);
for i=1:dim1
    for j=1:dim2
    for k=1:dim3
    SI=squeeze(img(i, j, k, :));
    SI=SI(1:2);
    if (max(SI)<SI_thres)
      tmp(i, j, k)=0.0;
    else
    [q w]=polyfit(b/1000, (log(SI))', 1);
    
% Use b unit: 1000s/mm2

% Negative ADC pixels are set to zero
        if q(1)<0.0
        tmp(i, j, k)=-1.*q(1);
        else
        tmp(i, j, k)=0.0;
        end
    end
    end
    end
end

WriteBrikEZ(img,info,'reorder_ADC',[prefix,'_ro']);

WriteBrikEZ(tmp,info,'reorder_ADC',[prefix,'_ADC']);


