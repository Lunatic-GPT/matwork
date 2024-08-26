function ress=mtimes(b,bb_all)

       

tinit=tic;

  

    load(b.fmap_ind);
    if b.adjoint
        ress = zeros([b.imSize,nmaps]);
    else
        ress = zeros([b.dataSize(1),Nc]);
    end
    
        a.dataSize=[b.dataSize(1),1];  %dont do a=b because you won't be able to access them in another function. 
        a.imSize=b.imSize;
        
        a.Nd = b.Nd;
        a.Jd = b.Jd;
        a.Kd = b.Kd;
        a.n_shift = b.n_shift;
        a.w=b.w;
        a.fmap_ind=b.fmap_ind;

        a.om = b.om(1:b.dataSize(1),:);
        
        if b.adjoint
            bb = bb_all( 1:b.dataSize(1),:);
        else
            bb=bb_all;
        end
        
        ress_tmp= NUFFT3D_worker(a,bb,b.maps,b.weights,b.adjoint);
        
        if b.adjoint
            ress=ress+ress_tmp;
        else
            ress(1:b.dataSize(1),:) = ress_tmp;
        end 
      


if b.adjoint
    disp('Total nufft3d mtimes adjoint time');
else
    disp('Total nufft3d mtimes time');
end
toc(tinit);


function ress=NUFFT3D_worker(a,bb,maps,weights,adjoint)


if adjoint
   nd=2;
else
    nd=4;
end
load(a.fmap_ind);  %ind1, ind2, and nro, Nc
    
if adjoint  % result is image
    
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
    
    
    clear maps weights ress;
    
    ress=ress2;
    toc;
   

else  %result is k-space data
    
%     tic;
%     disp('get maps from disk');
%     load(fmap_ind);  %ind1, ind2, and nro,Nc
%     [maps,weights]=get_maps_weights(prefix_map,a.imSize,nro,Nc,ind1,ind2,ncomp);
%     toc;
%     
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
% 
% disp('save worker output');
% tic;
% filename = [strtok(f_nufft,'.'),'_res.mat'];
% save(filename,'ress');
% 
% toc;
% 
% disp('Total time for worker');
% toc(tinit);



