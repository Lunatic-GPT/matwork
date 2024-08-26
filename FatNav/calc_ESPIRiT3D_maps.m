function calc_ESPIRiT3D_maps(kCalib_file,indList,nmaps,ine)
% kCalib_file should contain:
% imSize: nro*x*y
% 
% kCalib: npe*nro*npar*nch
% in the kspace for the 1st and 3rd dimensions.


load(kCalib_file);

Nc=size(kCalib,4);
kSize=[7,7];
eigThresh_im = 0.95;
eigThresh_k = 0.0085;
if size(kCalib,ine)~=imSize(ine)
    error('Size mismatch');
end

if ine==1
maps=zeros([length(indList),imSize(2),imSize(3),Nc,nmaps],'single');
weights=zeros([length(indList),imSize(2),imSize(3),nmaps],'single');
else
maps=zeros([imSize(1),imSize(2),length(indList),Nc,nmaps],'single');
weights=zeros([imSize(1),imSize(2),length(indList),nmaps],'single');   
end

troot=tic;
for i=1:length(indList)
    
    if ine==1
        kCalib_tmp=squeeze(kCalib(indList(i),:,:,:));
    else
        kCalib_tmp=squeeze(kCalib(:,:,indList(i),:));
    end
    [kernel,s] = dat2Kernel(kCalib_tmp,kSize);
    idx = max(find(s >= s(1)*eigThresh_k));
    if idx>200% most likely no signal; just noise
       % continue;
    end
  %  disp(indList(i));
  kernel_all(:,:,:,:,i)=kernel;
  if ine==1
    [M,W] = kernelEig(kernel(:,:,:,1:idx),imSize([2,3]));
    maps(i,:,:,:,:) = M(:,:,:,end-1:end);
    weights(i,:,:,:) = double((W(:,:,end-1:end) >  eigThresh_im));
  else
     [M,W] = kernelEig(kernel(:,:,:,1:idx),imSize([1,2]));
    maps(:,:,i,:,:) = M(:,:,:,end-1:end);
    weights(:,:,i,:) = double((W(:,:,end-1:end) >  eigThresh_im));  
  end
    
    fprintf('calc_ESPIRiT3D_maps: %d/%d - ind = %d; estimated remaining/total time: %d s/%d s\n',...
            i,length(indList),indList(i),round(toc(troot)/i*(length(indList)-i)),round(toc(troot)/i*length(indList)));
end

prefix=strtok(kCalib_file,'.');
name=sprintf('%s_maps_%d_%d.mat',prefix,indList(1),indList(end));
save(name, 'maps','weights','kernel_all','s','-v7.3');



