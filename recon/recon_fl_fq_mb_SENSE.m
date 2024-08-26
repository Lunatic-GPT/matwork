function recon_fl_fq_mb(fmb,fsb)

%fmb='meas_MID42_fl_fqZ_AF3_FID10324.dat';
%fsb='meas_MID44_fl_fqZ_SB_FID10326.dat';

% SENSE recon


[dsb,lin_sb,par_sb,sl_sb]=readMeasDat(fsb,inf,0,true);

[dmb,lin_mb,par_mb,sl_mb]=readMeasDat(fmb,inf,0,true);

prefix_mb=strtok(fmb,'.');
prefix_sb=strtok(fsb,'.');


nsl=readsPar([prefix_mb,'.pro'],'alFree[4]');
af=readsPar([prefix_mb,'.pro'],'alFree[6]');
nave=readsPar([prefix_mb,'.pro'],'lAverages');
nave_sb=readsPar([prefix_sb,'.pro'],'lAverages');
lin2=lin_mb(1:64*nave:end);
par2=par_mb(1:64*nave:end);

nro=size(dmb,1);
dmb=reshape(dmb,[nro,32,2,nave,length(dmb(:))/nro/64/nave]);
dsb=reshape(dsb,[nro,32,2,nave_sb,length(dsb(:))/nro/64/nsl/nave_sb,nsl]);

%%

dmb2=zeros(nro,32,2,max(lin2(:))+1,'single');

dsb2=zeros(nro,32,2,max(lin2(:))+1,nsl,'single');

    for i=1:length(lin2)
        dmb2(:,:,:,lin2(i)+1)=mean(dmb(:,:,:,:,i),4);
        dsb2(:,:,:,lin2(i)+1,:)=mean(dsb(:,:,:,:,i,:),4);
    end

dsb2=interleave2linear(dsb2,5);

nline=length(lin2);

fall=0*dsb2;
for i=1:nsl
    
f=exp(1i*(0:nline-1)*2*pi/af*(i-1-(nsl-1)/2));
f=reshape(f,[1,1,1,nline]);

fall(:,:,:,:,i)=repmat(f,[nro,32,2,1]);

end

%%

fdsb2=fft1c(dsb2,1);
fdmb2=fft1c(dmb2,1);


tmp=sos(fdmb2,2);
[tmp2,ind_max]=max(tmp(:));
ind_max=ind2subb(size(tmp),ind_max);
negphase=ind_max(1)-1;


im_sb=zeros(nro,nline,nsl,2,32,'single');
im_mb=zeros(nro,nline,2,32,'single');

        for k=1:32
            
            tmp=squeeze(fdsb2(:,k,:,:,:).*fall(:,k,:,:,:));
            tmp=permute(tmp,[3,1,4,2]);
            tmp=partialFT(tmp(:,1:negphase+nro/2,:,:),negphase);
     
            im_sb(:,:,:,:,k)=permute(tmp,[2,1,3,4]);
        %{  
            tmp=squeeze(fdmb2(:,k,:,:));
            tmp=permute(tmp,[3,1,2]);
            tmp=partialFT(tmp(:,1:negphase+nro/2,:),negphase);
            
            im_mb(:,:,:,k)=permute(tmp,[2,1,3]);
          %}  
            disp(k);
        end


im_recon=zeros(nro,nline,nsl,2,'single');

thr=0.08*max(vec(sos(im_sb,5)));

clust=clusterize2_2d(sos(im_sb(:,:,:,1,:),5)>thr,20);
clust2=clusterize2_2d_hole(clust,20);

gfactor=ones(nro,nline,nsl);

for i=1:nro
    for j=1:nline
        for k=1:2

            
            ref_k=1;
            ssel=squeeze(clust2(i,j,:));
            
            if ~any(ssel)
                continue;
            end
            
            norm=sos(im_sb(i,j,ssel,ref_k,:),5);
            
            mat=im_sb(i,j,ssel,ref_k,:)./repmat(norm,[1,1,1,1,32]);
        
            if sum(ssel)>1 
             mat=permute(squeeze(mat),[2,1]);
            else
                mat=squeeze(mat);
            end
         y=im_mb(i,j,k,:);
         
         im_recon(i,j,ssel,k)=mat\y(:);

         cov=mat'*mat;
         cov2=diag(diag(cov));
         
         gfactor(i,j,ssel)=sqrt(diag(inv(cov))./diag(inv(cov2)));
         
        end
    end
end

fall2=reshape(conj(fall(1,1,1,:,:)),[1,nline,nsl,1]);

fim=fft1c(im_recon,2).*repmat(fall2,[nro,1,1,2]);

im_mb2=ifft1c(fim,2);%.*repmat(clust2,[1,1,1,2]);

im_mb2=fftshift(fftshift(im_mb2,1),2);

%%
pc_mb=angle(im_mb2(:,:,:,1).*conj(im_mb2(:,:,:,2)))*180/pi;

%ph_sb=angle(im_sb2(:,:,:,1,:).*conj(im_sb2(:,:,:,2,:)))*180/pi;


save(['recon_',prefix_mb],'im_mb2');

save(['reconPC_',prefix_mb],'pc_mb');

save(['gfactor_',prefix_mb],'gfactor');

tmp=abs(im_mb2);
save(['reconMag_',prefix_mb],'tmp');

tmp=angle(im_mb2)*180/pi;
save(['reconPh_',prefix_mb],'tmp');

save(['recon_',prefix_sb],'im_sb');






