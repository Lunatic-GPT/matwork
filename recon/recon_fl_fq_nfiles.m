function recon_fl_fq_nfiles(fsb,interp_factor)

% if nphase is Inf; then set it equal to number of averages and do not use
% the physio data for retro-gating.

dsb2=[];
    for i=1:length(fsb)

        [dsb2_tmp,dsb1_32,pat,lin,npe]=sort_data(fsb{i});
        dsb2=cat(6,dsb2,dsb2_tmp);
        
    end
    
 
    dsb2=mean(dsb2,6);
%%



if pat>1
    dsb2_mean=mean(dsb2(:,:,:,:,1,:),6);
    lin2_unq=sort(unique(lin(:)));
    for i=1:size(dsb2,6)
        %dsb2(:,:,:,:,:,i)
        [im,dsb2(:,:,:,:,:,i)]=recon_grappa2D(fft1c(dsb2(:,lin2_unq+1,:,:,:,i),1),lin2_unq+1,npe,fft1c(dsb2_mean(:,lin2_unq+1,:,:),1),lin2_unq+1); 
    end
else
    
    dsb2= fft1c(dsb2,1); 
end
%%

tmp=sos(dsb2(:,:,:,:,1,1),4);
[tmp2,ind_max]=max(tmp(:));
ind_max=ind2subb(size(tmp),ind_max);
negphase=ind_max(1)-1;
%%
sz=size(dsb2);
im_sb=single(zeros(sz));

for i=1:size(dsb2,4)
    for j=1:size(dsb2,5)
        for k=1:size(dsb2,3)
            for iphase=1:size(dsb2,6)
                tmp=dsb2(:,:,k,i,j,iphase);
               
               if abs((sz(1)/2-negphase)/sz(1))>0.1 % no partial Fourier recon if less than 10%
                    tmp=cat(1,zeros(sz(1)/2-negphase,size(tmp,2)),tmp(1:negphase+sz(1)/2,:));
                    im_sb(:,:,k,i,j,iphase)=partialFT_pocs(tmp,40,true);

             %       tmp=partialFT_pocs(tmp,40,true);
              %      tmp=fft1c(fft1c(tmp,1),2);
               %     im_sb(:,:,k,i,j,iphase)=ifft1c(ifft1c(tmp,1),2);       
                else 
                    im_sb(:,:,k,i,j,iphase)=ifft1c(ifft1c(tmp,1),2);
                end 
            end 
        end        
    end
end

%%


%mkdir(prefix);
mid='';
for i=1:length(fsb)
    mid_tmp=strtok_no(fsb{i},'_',2);
    mid=[mid,mid_tmp,'_'];
end

mid=mid(1:end-1);
%dsb1_32=dsb(:,1:32);
% save(['recon_',mid,'.mat'],'im_sb','dsb1_32','freePara');
 
do_interp_dim12(im_sb,interp_factor,mid,dsb1_32,['recon_',mid]);



function [dsb2,dsb1_32,pat,lin,npe]=sort_data(fsb)

prefix=strtok(fsb,'.');
prefix_dir=fullfile(prefix,prefix);
if ~exist([prefix,'.mat'],'file')
[dsb,lin,par,sl,ushSet,timeStamp,freePara]=readMeasDat(fsb,inf,0,true);
save([prefix,'.mat'],'dsb','lin','par','sl','ushSet','freePara');
%nois=dsb(:,1:32);
%save([prefix_dir,'.mat'],'nois');  % noise data for later use.

else
    
load([prefix,'.mat']);

if exist('Set','var')  % new format
ushSet=Set;
lin=Line;
dsb=Data;
clear Data;
end



end

dsb1_32=dsb(:,1:32);

nave=readsPar([prefix_dir,'.pro'],'lAverages');

seg=readsPar([prefix_dir,'.pro'],'lSegments');


if isempty(nave)
    nave=1;
end

pat=readsPar([prefix_dir,'.pro'],'lAccelFactPE');
npe=readsPar([prefix_dir,'.pro'],'lPhaseEncodingLines');
%%

%tr=readsPar([prefix_dir,'.pro'],'alTR[0]');


%for compatibility with older files


nvenc=max(ushSet)+1;

%%


  dsb2=reorder_fl_fq_data(dsb,lin,prefix_dir,nvenc,freePara(:,4),Inf,[],[],false,90,[],[]);
  






  
  
