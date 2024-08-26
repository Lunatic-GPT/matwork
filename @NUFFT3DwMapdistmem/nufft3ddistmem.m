function  res = nufft3ddistmem(k,w,shift,imSize,nblock,prefix_map,fmap_ind,ine,distMem)
%
if ~exist('nblock','var')
    nblock=1;
end
    res.distMem=distMem;
    res.Nd = imSize;
    res.Jd = [5,5,5];
    res.Kd = floor([imSize*1.4]);
    res.n_shift = imSize/2 + shift;
  %  res.st = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');  calculate on
  %  fly to save memory.
    res.nblock=nblock;
    res.om=k*2*pi;
    res.adjoint = 0;
   
    res.imSize = imSize;
    res.dataSize = [size(k,1),1];
    res.w = sqrt(w);
    res.prefix_map = prefix_map;
    res.fmap_ind = fmap_ind;
    
    if ~distMem 
    disp('get maps from disk');
    load(fmap_ind);
    [res.maps,res.weights]=get_maps_weights(prefix_map,imSize,nro,Nc,ind1,ind2,nmaps,ine);
    disp('get map time');
    end
    
    res = class(res,'nufft3ddistmem');

    
function [maps,weights]=get_maps_weights(prefix,imSize,nro,Nc,ind1,ind2,nmaps,ine)


maps=zeros([nro,imSize(2),imSize(3),Nc,nmaps],'single');
weights=ones([nro,imSize(2),imSize(3),nmaps],'single');

for i=1:length(ind1)
    while 1
        try
            tmp=load(sprintf('%s_%d_%d.mat',prefix,ind1(i),ind2(i)));
            break;
        catch
            fprintf('Reading %s_%d_%d.mat error; try again!\n',prefix,ind1(i),ind2(i));
            continue;
        end
    end
    if ine==1
    maps(ind1(i):ind2(i),:,:,:,:)=tmp.maps(:,:,:,:,end-nmaps+1:end);
    weights(ind1(i):ind2(i),:,:,:)=tmp.weights(:,:,:,end-nmaps+1:end);
    else
      maps(:,:,ind1(i):ind2(i),:,:)=tmp.maps(:,:,:,:,end-nmaps+1:end);
      weights(:,:,ind1(i):ind2(i),:)=tmp.weights(:,:,:,end-nmaps+1:end); 
    end
end

