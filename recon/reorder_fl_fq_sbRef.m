function dsb2=reorder_fl_fq_sbRef(dsb,lin,prefix,nvenc,nsl,order_slice,Slice)

% the output is of dimension [nro,lPhaseEncodingLines,nsl,32,nvenc,nphase];
if ~exist('order_slice','var')
    order_slice=true;
end

nvenc=double(nvenc);
nave=readsPar([prefix,'.pro'],'lAverages');


seg=readsPar([prefix,'.pro'],'lSegments');

if isempty(nave)
    nave=1;
end

pat=readsPar([prefix,'.pro'],'lAccelFactPE');
nsl_first=readsPar([prefix,'.pro'],'lConc');
npe=readsPar([prefix,'.pro'],'lPhaseEncodingLines');
if pat>1
    lin=lin(33:end);
    dsb=dsb(:,33:end);
end

 
    lin2=lin(1:32:end);
    
    nphase=1;
    
%lin2=lin2(1:end/nsl);

  ivenc = repmat(1:nvenc,[1,ceil(length(lin2)/nvenc)]);
  ivenc=ivenc(1:length(lin2));

    
  if ~exist('Slice','var')
      Slice = repmat(1:nsl,[length(lin2)/nsl,1]);
      Slice=Slice(:);
  else
  Slice = Slice+1;
  Slice=Slice(33:32:end);
  end
  
%lin2=reshape(lin2,[seg,nave,length(lin2)/seg/nave]);
%lin2=reshape(lin2(:,1,:),[seg,length(lin2(:))/seg/nave]);



nro=size(dsb,1);

%%

    dsb2=zeros(nro,32,nvenc,npe,nsl,'single');
    
  %  dsb=reshape(dsb,[nro,32,nvenc,seg,nave,length(dsb(:))/nro/32/nvenc/nsl/nave/seg,nsl]);
   
        dsb=reshape(dsb,[nro,32,size(dsb,2)/32]);
   
    lin2_unq=unique(lin2);
    
    for k=1:nsl
    for i=1:length(unique(lin2))
     for j=1:nvenc
        dsb2(:,:,j,lin2_unq(i)+1,k)=mean(dsb(:,:,lin2==lin2_unq(i) & ivenc==j&Slice == k),3);
     end
    
    end
    
    end
  %{  
    for j=1:seg
        for i=1:size(lin2,2)
            if pat==1
            dsb2(:,:,:,lin2(j,i)+1,:)=mean(dsb(:,:,:,j,:,i,:),5);
            else
            dsb2(:,:,:,lin2(j,i)+1,:)=mean(dsb(:,:,:,j,:,i,:),5);
            
            end
        end
    end
%}
dsb2=permute(dsb2,[1,4,5,2,3,6]);

if nsl==1 || nsl_first <nsl % the reference acquired in multiple scans. do not reorder slices in this case
    return;
end


if order_slice
    
    ucMode = readsPar([prefix,'.pro'],'sSliceArray.ucMode');
    
    if strcmp(ucMode,'0x4')  % interleaved
        dsb2=interleave2linear(dsb2,3);
    elseif strcmp(ucMode,'0x2')  %descend
        dsb2=flip(dsb2,3);
    end
     
end





  
  
