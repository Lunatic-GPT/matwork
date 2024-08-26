function dsb2=reorder_fl_fq_data_phaseStabOn(dsb,lin,prefix,nvenc,physio,nphase,do_pc)

% the output is of dimension [nro,lPhaseEncodingLines,nsl,32,nvenc,nphase];
% only tested for pat = 1 so far
mid=strtok_no(prefix,'_',2);

nvenc=double(nvenc);
nave=readsPar([prefix,'.pro'],'lAverages');

seg=readsPar([prefix,'.pro'],'lSegments');

if isempty(nave)
    nave=1;
end

pat=readsPar([prefix,'.pro'],'lAccelFactPE');
nsl=readsPar([prefix,'.pro'],'lConc');
tr=readsPar([prefix,'.pro'],'alTR[0]');
npe=readsPar([prefix,'.pro'],'lPhaseEncodingLines');

try
    
nsl_mb=readsPar([prefix,'.pro'],'alFree[24]');  %% was alFree[5]
af=readsPar([prefix,'.pro'],'alFree[6]');

catch
    nsl_mb=1;
    af=1;
end

if nsl_mb==1
    af=1;
end

nro=size(dsb,1);

%if phaseStab
    lin=lin(32*nvenc*2+1:end);
    
    lin=lin(1:2:end);
    
    physio=physio(32*nvenc*2+1:end);
    
    physio=physio(1:2:end);
    
    dsb=dsb(:,32*nvenc*2+1:end);
    
    
    dsb=reshape(dsb,[size(dsb,1),32,2,size(dsb,2)/64]);
    
    dsb_ref=dsb(:,:,2,:);
    dsb=dsb(:,:,1,:);
    
if exist('do_pc','var') && do_pc 
    
    dsb_ref=reshape(dsb_ref,[size(dsb_ref,1),size(dsb_ref,2),size(dsb_ref,3),nave*nvenc,nsl_mb*af,size(dsb_ref,4)/nave/nvenc/nsl_mb/af]);
    
dsb_ref_mn=mean(mean(dsb_ref,4),6);

dsb_ref=mean(mean(dsb_ref.*conj(repmat(dsb_ref_mn,[1,1,1,size(dsb_ref,4),1,size(dsb_ref,6)])),2),1);
dsb_ref=dsb_ref./abs(dsb_ref);
dsb_ref=reshape(dsb_ref,[1,1,1,length(dsb_ref(:))]);

dsb=dsb.*repmat(conj(dsb_ref),[nro,32,1,1]);

end

%end

% 
% if pat>1
%     lin=lin(33:end);
%     dsb=dsb(:,33:end);
%     physio=physio(33:end);
% end

if nphase==1
  
    lin2=lin(1:32*nvenc:end);  
  
    lin2=lin2(1:end/nsl);
    lin2=reshape(lin2,[seg,nave,length(lin2)/seg/nave]);
    lin2=reshape(lin2(:,1,:),[seg,length(lin2(:))/seg/nave]);
else

    lin2=lin(1:32:end);
    physio=physio(1:32:end);

    lgrad=0*lin2;
    
    for i=1:nvenc
     lgrad(i:nvenc:end)=i;
    end
    
    %%divide into phases based on peaks only
    %ph=physioTrigger(physio,0.43,tr*1e-6/seg/nvenc,nphase,0.1,0,1);

    %% divide into phases based both on physio peaks and k-space data
    
    if exist(['FakePeaks_',mid,'.txt'],'file')
        
        fakePeaks=load(['FakePeaks_',mid,'.txt']);
    else
        fakePeaks=[];
    end
    
    ph=physioTrigger(physio,0.43,tr*1e-6/seg/nvenc,nphase,0.1,0,1,fakePeaks);
    
 %   pause;
    % check whether any k-space line not acquired
    
    lin2_unq=unique(lin2);
    nacq=zeros(nvenc,length(lin2_unq),nphase);
    
   for i=1:length(lin2_unq)
       for j=1:nphase
        for k=1:nvenc
         nacq(k,i,j) = sum(lin2==lin2_unq(i)&lgrad==k&ph==j);        
        end
       end
   end  
    
   fprintf('Empty lines (%f%%):\n ',100*sum(nacq(:)==0)/length(nacq(:)));
   ind_tmp=ind2subb(size(nacq),find(nacq==0));
   fprintf('venc   line     phase\n');
   disp(ind_tmp);
         
end


%%

    dsb2=zeros(nro,32,nvenc,npe,nsl,nphase,'single');
if ~do_retro
    
    dsb=reshape(dsb,[nro,32,nvenc,seg,nave,length(dsb(:))/nro/32/nvenc/nsl/nave/seg,nsl]);
   
    
    for j=1:seg
        for i=1:size(lin2,2)
            if pat==1
            dsb2(:,:,:,lin2(j,i)+1,:)=mean(dsb(:,:,:,j,:,i,:),5);
            else
            dsb2(:,:,:,lin2(j,i)+1,:)=mean(dsb(:,:,:,j,:,i,:),5);
            
            end
        end
    end
else
    if nsl>1
        error('Multiple slices not possible with retro-gating');
    end
    
    dsb=reshape(dsb,[nro,32,length(dsb(:))/32/nro]);
    
    for i=1:length(lin2_unq)
        for j=1:nphase
            for k=1:nvenc
                if nacq(k,i,j)>0
                     dsb2(:,:,k,lin2_unq(i)+1,1,j) = mean(dsb(:,:,lin2==lin2_unq(i)&lgrad==k&ph==j),3);
                   
                        
                else
                    ind_tmp= (lin2==lin2_unq(i)&lgrad==k&ph==j-1) | (lin2==lin2_unq(i)&lgrad==k&ph==j+1);
                    
                    if sum(ind_tmp)>0
                        dsb2(:,:,k,lin2_unq(i)+1,1,j)=mean(dsb(:,:,ind_tmp),3);
                    else
                        ind_tmp= (lin2==lin2_unq(i)&lgrad==k&ph==j-2) | (lin2==lin2_unq(i)&lgrad==k&ph==j+2);
                        
                        if sum(ind_tmp)>0
                            dsb2(:,:,k,lin2_unq(i)+1,1,j)=mean(dsb(:,:,ind_tmp),3);
                            warning('k-space line moment %d phase %d line %d set to 2nd nearest phase',k,j,i);
                        else
                            warning('k-space line moment %d phase %d line %d set to zero',k,j,i);
                        end
                    end
                end
            end
        end
    end
    
end


dsb2=permute(dsb2,[1,4,5,2,3,6]);

if nsl==1
    return;
end


ucMode = readsPar([prefix,'.pro'],'sSliceArray.ucMode');

if strcmp(ucMode,'0x4')  % interleaved
    dsb2=interleave2linear(dsb2,3);
elseif strcmp(ucMode,'0x2')
    dsb2=flip(dsb2,3);
end










  
  
