function [dsb2,hr]=reorder_fl_fq_data(dsb,lin,pro,nvenc,physio,nphase,physio_a,plot_results,maxRate)

% reorder_fl_fq_data(dsb,Line,prefix_dir,nvenc,freePara(:,4),nphase,f_fake,f_MON,true,maxRate,f_real,f_interp);
% the output is of dimension [nro,lPhaseEncodingLines,nsl,32,nvenc,nphase];
% 7/25/2017: dsb2=reorder_fl_fq_data(dsb,lin,prefix,nvenc,do_retro,physio,nphase)
% changed to dsb2=reorder_fl_fq_data(dsb,lin,prefix,nvenc,physio,nphase,f_fake,f_MON)
% 10/13/2020: input argument changed from prefix to pro

if isfield(physio_a,'f_fake')
    f_fake=physio_a.f_fake;
else
    f_fake=[];
end

if isfield(physio_a,'f_MON')
    f_MON=physio_a.f_MON;
else
    f_MON=[];
end

if isfield(physio_a,'f_real')
    f_real=physio_a.f_real;
else
    f_real=[];
end

if isfield(physio_a,'f_interp')
    f_interp=physio_a.f_interp;
else
    f_interp=[];
end

nvenc=double(nvenc);
nave=readsPar(pro,'lAverages');

seg=readsPar(pro,'lSegments');

if isempty(nave)
    nave=1;
end



pat=readsPar(pro,'lAccelFactPE');
nsl=readsPar(pro,'lConc');
tr=readsPar(pro,'alTR[0]');
npe=readsPar(pro,'lPhaseEncodingLines');

 tr=tr/1e6;
if pat>1
    lin=lin(33:end);  % the noise scan is before the first FatNav scan
    dsb=dsb(:,33:end);
    physio=physio(33:end);
end

if length(lin)<length(physio)
   fatnav=true;
else
    fatnav=false;
end


if fatnav
    tr_nav=readsPar(pro,'sWiPMemBlock.adFree[3]');
    trimg_nav=readsPar(pro,'adFree[5]');
    dummy=readsPar(pro,'alFree[33]');
    trInterval=readsPar(pro,'alFree[30]');
    ncal = 984;
    nnav=384;
    tr_nav=tr_nav/1e3;
    trimg_nav=trimg_nav/1e3;
   
    t=[];
    indt=[];
    for i=0:length(unique(lin))-1
        
       if mod(i,trInterval)==0
           
           if i==0
             t0=0;
             t(end+1:end+ncal) = t0+(1:ncal)*tr_nav;
             t0=t0+(ncal-nnav)*tr_nav;
           else
               t0=t(end);
               t(end+1:end+nnav) = t0+(1:nnav)*tr_nav;
           end
           t(end+1:end+dummy*nvenc)=t0+trimg_nav+(1:dummy*nvenc)*tr/nvenc;
       end
       
          t(end+1:end+nave*nvenc)=t(end)+tr/nvenc*(1:nave*nvenc);
          indt(end+1:end+nave*nvenc)=length(t)-nave*nvenc+1:length(t);
    
    end
end



if nphase==1 || isinf(nphase)
  
    lin2=lin(1:32:end);  

    lin2=lin2(1:end/nsl);
    
    iseg=repmat(1:seg,[nvenc,1]);
    
        
    iseg=repmat(iseg(:),[ceil(length(lin2)/length(iseg(:))),1]);
    
    iseg=iseg(1:length(lin2));
    
    ivenc = repmat(1:nvenc,[1,ceil(length(lin2)/nvenc)]);
    ivenc=ivenc(1:length(lin2));
    
    
%    lin2=reshape(lin2,[seg,nave,length(lin2)/seg/nave]);
%    lin2=reshape(lin2(:,1,:),[seg,length(lin2(:))/seg/nave]);
    

    
else

    lin2=lin(1:32:end);
    physio=physio(1:32:end);

    lgrad=0*lin2;
    
    for i=1:nvenc
     lgrad(i:nvenc:end)=i;
    end
    
    
    fakePeaks =load_text_number(f_fake);
    
    Peaks2MON = load_text_number(f_MON);
    
    realPeaks = load_text_number(f_real);
    
    interpPeaks = load_text_number(f_interp);
    
minInterval =60/maxRate;
if ~fatnav
    [ph,tmp,hr]=physioTrigger(physio,minInterval,tr/2,nphase,0.1,0,plot_results,fakePeaks,Peaks2MON,realPeaks);  
 %  [ph,trig]=physioTrigger(data,minInterval,TR,nphase,trig_delay,thr,plot_results,fakePeaks,Peaks2MON)
else
    
 indt=indt(1:length(lin2));
   
 t=t(1:length(physio));%in case the end of the data file is corrupt;
    ph=physioTriggerNew(physio,minInterval,t,nphase,0.1,0,plot_results,fakePeaks,Peaks2MON,realPeaks,interpPeaks);
    
ph=ph(indt);
end
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
    
%    if sum(nacq(:)==0)>0
%    fprintf('Empty lines (%f%%):\n ',100*sum(nacq(:)==0)/length(nacq(:)));
%    ind_tmp=ind2subb(size(nacq),find(nacq==0));
%    fprintf('venc   line     phase\n');
%    disp(ind_tmp);
%          
%    else
%        disp('Nice! No un-filled k-space lines!');
%    end
% 


end
%%

nro=size(dsb,1);

if isinf(nphase)
   
    nphase_real=sum(lin2==lin2(1))/nvenc;
    dsb2=zeros(nro,32,nvenc,npe,nsl,nphase_real,'single');

else
    dsb2=zeros(nro,32,nvenc,npe,nsl,nphase,'single');
end

    totalLine=ceil(size(dsb,2)/32/nsl)*32*nsl;
    
  if size(dsb,2)-totalLine<0
     warning('%d lines missing in the data; Fill with zero',totalLine-size(dsb,2));  
         dsb(:,end+1:totalLine)=0;  % fill the missing lines with 0;
  end

if nphase==1 || isinf(nphase)
    
   % dsb=reshape(dsb,[nro,32,nvenc,seg,nave,length(dsb(:))/nro/32/nvenc/nsl/nave/seg,nsl]);
    dsb=reshape(dsb,[nro,32,size(dsb,2)/32/nsl,nsl]);
    
    lin2_unq=unique(lin2);
    
    for i=1:length(unique(lin2))
     for j=1:nvenc
        if nphase==1
           dsb2(:,:,j,lin2_unq(i)+1,:,:)=mean(dsb(:,:,lin2==lin2_unq(i) & ivenc==j,:),3);
        else
            sel=lin2==lin2_unq(i) & ivenc==j;
            nsel=sum(sel);
            if nsel~=size(dsb2,6)
               warning('Mismatch for line %d venc %d: expect = %d, got = %d. Fill zero.',lin2_unq(i),j,size(dsb2,6),sum(sel)); 
            end
           dsb2(:,:,j,lin2_unq(i)+1,:,1:nsel)=permute(dsb(:,:,sel,:),[1,2,4,3]);
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


ucMode = readsPar(pro,'sSliceArray.ucMode');

if strcmp(ucMode,'0x4')  % interleaved
    dsb2=interleave2linear(dsb2,3);
elseif strcmp(ucMode,'0x2')
    dsb2=flip(dsb2,3);
end










  
  
