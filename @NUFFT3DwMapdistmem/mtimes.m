function ress=mtimes(b,bb_all)

       

tinit=tic;
if b.distMem
    ress=mtimes_distMem(b,bb_all);
else
    
    nblock=b.nblock;
    ind=linspace(0,b.dataSize(1),nblock+1);
    ind=round(ind);
%     
%     tic;
%     disp('get maps from disk');
    load(b.fmap_ind);
%     [maps,weights]=get_maps_weights(b.prefix_map,b.imSize,nro,Nc,ind1,ind2,nmaps);
%     disp('get map time');
%     toc;
    
    if b.adjoint
        ress = zeros([b.imSize,nmaps]);
    else
        ress = zeros([b.dataSize(1),Nc]);
    end
    
    for i=1:nblock
        a.dataSize=[ind(i+1)-ind(i),1];  %dont do a=b because you won't be able to access them in another function. 
        a.imSize=b.imSize;
        
        a.Nd = b.Nd;
        a.Jd = b.Jd;
        a.Kd = b.Kd;
        a.n_shift = b.n_shift;
        a.w=b.w;
        a.fmap_ind=b.fmap_ind;

        a.om = b.om(ind(i)+1:ind(i+1),:);
        
        if b.adjoint
            bb = bb_all( ind(i)+1:ind(i+1),:);
        else
            bb=bb_all;
        end
        
        ress_tmp= NUFFT3DNonDistMem_worker(a,bb,b.maps,b.weights,b.adjoint);
        
        if b.adjoint
            ress=ress+ress_tmp;
        else
            ress(ind(i)+1:ind(i+1),:) = ress_tmp;
        end 
    end  
end

if b.adjoint
    disp('Total nufft3d mtimes adjoint time');
else
    disp('Total nufft3d mtimes time');
end
toc(tinit);
% 
% function [maps,weights]=get_maps_weights(prefix,imSize,nro,Nc,ind1,ind2,nmaps)
% 
% maps=zeros([imSize(1),nro,imSize(3),Nc,nmaps],'single');
% weights=ones([imSize(1),nro,imSize(3),nmaps],'single');
% 
% for i=1:length(ind1)
%     tmp=load(sprintf('%s_%d_%d.mat',prefix,ind1(i),ind2(i)));
%     maps(:,ind1(i):ind2(i),:,:,:)=tmp.maps(:,:,:,:,end-nmaps+1:end);
%     weights(:,ind1(i):ind2(i),:,:)=tmp.weights(:,:,:,end-nmaps+1:end);
% end


function ress = mtimes_distMem(b,bb_all)
% performs the normal nufft

tinit=tic;



nblock=b.nblock;
ind=linspace(0,b.dataSize(1),nblock+1);
ind=round(ind);
fres={};
fid=fopen('bsub_NUFFT3DDistMem_worker.sh','w');
for i=1:nblock
    
    tic;
    f_nufft= sprintf('nufft_init_%dblocks_%d.mat',nblock,i);
    
    
    %  if ~exist(f_nufft,'file')
    disp(['saving ',f_nufft]);
    a.dataSize=[ind(i+1)-ind(i),1];
    a.imSize=b.imSize;
    a.om = b.om(ind(i)+1:ind(i+1),:);
    a.Nd = b.Nd;
    a.Jd = b.Jd;
    a.Kd = b.Kd;
    a.n_shift = b.n_shift;
    a.w=b.w;
    save(f_nufft,'a');
    toc;
    % end
    tic;
    if b.adjoint
        
        f_data = sprintf('nufft_3DDistMem_data_%dblocks_%d.mat',nblock,i);
        
        disp(['Save ', f_data]);
        bb = bb_all( ind(i)+1:ind(i+1),:);
        save(f_data,'bb');
        
    else
        f_data = sprintf('nufft_3DDistMem_data_%dblocks.mat',nblock);
        if i==1
            
            disp(['Save ', f_data]);
            bb=bb_all;
            save(f_data,'bb');
        end
        
    end
    
    toc;
    
    fres{i}= [strtok(f_nufft,'.'),'_res.mat'];
    
    if exist(fres{i},'file')
        delete(fres{i});
    end
    fprintf(fid,'bsub -M 40 -o lognufft3DDistMem%d.%%J matbgk "NUFFT3DDistMem_worker(''%s'',''%s'',''%s'',''%s'',%d)" lognufft3DDistMem_mat%d\n',i,b.prefix_map,f_nufft,f_data,b.fmap_ind,b.adjoint,i);
end

fclose(fid);

cmd='source bsub_NUFFT3DDistMem_worker.sh';
unix(cmd);

debug=false;
if ~debug
    tic;
    disp('waiting for files');
    wait_for_files(fres);
    toc;
    
    tic;
    disp('assemble the results');
    
end

load(b.fmap_ind);
if b.adjoint
    ress = zeros([b.imSize,nmaps],'single');
else
    ress = zeros([b.dataSize(1),Nc],'single');
end
if ~debug
    for i=1:length(fres)
        tmp = load(fres{i});
        if b.adjoint
            ress=ress+tmp.ress;
        else
            ress(ind(i)+1:ind(i+1),:) = tmp.ress;
        end
    end
end

ress=double(ress);% required otherwise, error for single divide by func which is sparse;
toc;
if b.adjoint
    disp('Total nufft3d mtimes adjoint time');
else
    disp('Total nufft3d mtimes time');
end
toc(tinit);

function wait_for_files(flist,nsecond)
found =false;

if ~exist('nsecond','var')
    nsecond = 60;
end

while ~found
    pause(nsecond);
    found=true;
    for i=1:length(flist)
        dir_str=dir(flist{i});
        if isempty(dir_str) || dir_str.bytes<1000
            found = false;
            fprintf('File %s not found\n',flist{i});
            break;
        else
            
        end
    end
    
end

