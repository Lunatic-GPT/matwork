function NUFFT3DDistMem_worker(prefix_map,f_nufft,f_data,fmap_ind,adjoint)


tinit=tic;

tic;
disp('load f_nufft');
load(f_nufft);% load structure "a" that contains nufft parameters
toc;

tic;
disp('load data');
load(f_data);% bb: kdata or imdata data
toc;


if adjoint
   nd=2;
else
    nd=4;
   
end




if adjoint
    
   
      load(fmap_ind);  %ind1, ind2, and nro, Nc
    
    disp('do nufft adjoint');
    
    if length(a.w(:))>1
        bb=bb.*repmat(a.w,[1,size(bb,nd)]);
    else
        bb=bb.*a.w(:,:);
    end
    

    ress = zeros([a.imSize,size(bb,nd)],'single');
    

   %%{ 
tic;
disp('nufft_init');
st = nufft_init(a.om, a.Nd, a.Jd, a.Kd, a.n_shift,'kaiser');
toc;

    for m=1:Nc
        tic;
        b = bb(:,m);

        res = nufft_adj(double(b), st)/sqrt(prod(a.imSize));
        ress(:,:,:,m)= reshape(res, a.imSize);
        fprintf('remaining = %d s; nufft_adj per chan = %f s\n',...
            round(toc*(Nc-m)),toc);
    end
    
    clear st res bb b;
  %}
       tic;
    disp('get maps from disk');
  
    [maps,weights]=get_maps_weights(prefix_map,a.imSize,nro,Nc,ind1,ind2,ncomp);
    toc;
    
    
    disp('do ESP adjoint');
    tic;
     weights=sqrt(weights);
     nv = size(weights,4);
     ress2 = zeros([a.imSize,nv],'single');
     
    for n=1:nv
        for j=1:size(maps,4)
          ress2(:,:,:,n) = ress2(:,:,:,n)+conj(maps(:,:,:,j,n)).*ress(:,:,:,j).*weights(:,:,:,n);
        end
    end
    
%     for j=1:size(ress,4)
%         
%         restmp = ress(:,:,:,j);
%         save(sprintf('tempChan_%d.mat',j),restmp);
%     end
%     
    
    clear maps weights ress;
    
    ress=ress2;
    toc;
   

else
    
    tic;
    disp('get maps from disk');
    load(fmap_ind);  %ind1, ind2, and nro,Nc
    [maps,weights]=get_maps_weights(prefix_map,a.imSize,nro,Nc,ind1,ind2,ncomp);
    toc;
    
    weights=sqrt(weights);
     
    nv = size(weights,4);
    tic;
    b = zeros([a.imSize,Nc],'single');
   for n=1:Nc
       for j=1:nv
       b(:,:,:,n) = b(:,:,:,n)+ maps(:,:,:,n,j).*weights(:,:,:,j).*bb(:,:,:,j);
       end
   end
    toc;
    
    clear maps weights bb;
     ress = zeros([size(a.om,1),Nc],'single');
%%{
     tic;

disp('nufft_init');
st = nufft_init(a.om, a.Nd, a.Jd, a.Kd, a.n_shift,'kaiser');
toc;
    
    for m=1:Nc
        tic;
   
        res = nufft(double(b(:,:,:,m)), st)/sqrt(prod(a.imSize));
        ress(:,m) = res(:);
        
        fprintf('remaining = %d s; nufft per chan = %f s\n',...
            round(toc*(Nc-m)),toc);
    end
    
    
    clear st res;
    %}
    clear b;
    %         end
    if length(a.w(:))>1
        ress=ress.*repmat(a.w,[1,size(ress,2)]);
    else
        
        ress=ress*a.w;
    end
    
    
end

disp('save worker output');
tic;
filename = [strtok(f_nufft,'.'),'_res.mat'];
save(filename,'ress');

toc;

disp('Total time for worker');
toc(tinit);

function [maps,weights]=get_maps_weights(prefix,imSize,nro,Nc,ind1,ind2,ncomp)

maps=zeros([imSize(1),nro,imSize(3),Nc,ncomp],'single');
weights=ones([imSize(1),nro,imSize(3),ncomp],'single');

for i=1:length(ind1)
    tmp=load(sprintf('%s_%d_%d.mat',prefix,ind1(i),ind2(i)));
    maps(:,ind1(i):ind2(i),:,:,:)=tmp.maps(:,:,:,:,end-ncomp+1:end);
    weights(:,ind1(i):ind2(i),:,:)=tmp.weights(:,:,:,end-ncomp+1:end);
end


