function dsb2=reorder_fl_fq_retro_data(dsb,lin,set,prefix,nvenc,nphase)

% the output is of dimension [nro,lPhaseEncodingLines,nsl,32,nvenc,nphase];

nvenc=double(nvenc);
nave=readsPar([prefix,'.pro'],'lAverages');

seg=readsPar([prefix,'.pro'],'lSegments');

if isempty(nave)
    nave=1;
end

pat=readsPar([prefix,'.pro'],'lAccelFactPE');
nsl=readsPar([prefix,'.pro'],'lConc');

npe=readsPar([prefix,'.pro'],'lPhaseEncodingLines');
if pat>1
    lin=lin(33:end);
    dsb=dsb(:,33:end);
    set=set(:,33:end);
end



    lin2=lin(1:32:end);
    set=set(1:32:end);

    
    
    dlin2=diff(lin2);
    
    ph=divide2Phases(find(dlin2~=0)+1,length(lin2),nphase,0.3);
    
    
   %  ph = divide2Phases(trig,npts,nphase,trig_delay)
     
 %   pause;
    % check whether any k-space line not acquired
    
    lin2_unq=unique(lin2);
    nacq=zeros(nvenc,length(lin2_unq),nphase);
    
   for i=1:length(lin2_unq)
       for j=1:nphase
        for k=1:nvenc
         nacq(k,i,j) = sum(lin2==lin2_unq(i)&set==k&ph==j);        
        end
       end
   end  
    
   fprintf('Empty lines (%f%%):\n ',100*sum(nacq(:)==0)/length(nacq(:)));
   ind_tmp=ind2subb(size(nacq),find(nacq==0));
   fprintf('venc   line     phase\n');
   disp(ind_tmp);
         



nro=size(dsb,1);

%%

    dsb2=zeros(nro,32,nvenc,npe,nsl,nphase,'single');
    if nsl>1
        error('Multiple slices not possible with retro-gating');
    end
    
    dsb=reshape(dsb,[nro,32,length(dsb(:))/32/nro]);
    
    for i=1:length(lin2_unq)
        for j=1:nphase
            for k=1:nvenc
                if nacq(k,i,j)>0
                     dsb2(:,:,k,lin2_unq(i)+1,1,j) = mean(dsb(:,:,lin2==lin2_unq(i)&set==k&ph==j),3);
                   
                        
                else
                    ind_tmp= (lin2==lin2_unq(i)&set==k&ph==j-1) | (lin2==lin2_unq(i)&set==k&ph==j+1);
                    
                    if sum(ind_tmp)>0
                        dsb2(:,:,k,lin2_unq(i)+1,1,j)=mean(dsb(:,:,ind_tmp),3);
                    else
                        ind_tmp= (lin2==lin2_unq(i)&set==k&ph==j-2) | (lin2==lin2_unq(i)&set==k&ph==j+2);
                        
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






function ph = divide2Phases(trig,npts,nphase,trig_delay)
 
% trig: indices for peaks; should be between 1 and npts 
% npts: number of physio signal time points
% nphase: number of phases to divide each cardiac cycle into

ph=zeros(1,npts);
trig=trig+trig_delay;
for i=1:npts
 
    itr=find(i>=trig(1:end-1)&i<trig(2:end));
    
  
    if isempty(itr)
        if i<trig(1) %before first peak
      
            trig0=2*trig(1)-trig(2);
            
          ph(i)=round((i-trig0)/(trig(2)-trig(1))*nphase);
          
        else  %after last peak
        
           if round((npts-trig(end))/(trig(end)-trig(end-1))*nphase)>nphase % too long
               ph(i)=round((i-trig(end))/(npts-trig(end))*nphase);
           else
             ph(i)=round((i-trig(end))/(trig(end)-trig(end-1))*nphase);
           end
        end
    else
      
        ph(i)=round((i-trig(itr))/(trig(itr+1)-trig(itr))*nphase);
        
    end
    
end

ph(ph==0)=nphase;
for i=1:nphase
    fprintf('Phase %d = %d\n',i,sum(ph==i));
end




  
  
