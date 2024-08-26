function d1bk=epi_phase_corr(d1b,ref)
% d1b k-space data; nro*npe*...
% ref: reference data; nro*2*...; 
    d1b=ifft1c(d1b,1);
    ref=ifft1c(ref,1);
    
    
    %d1b=ifft(d1b,[],1);
    %ref=ifft(ref,[],1);
    
    ph=angle(ref);
 
    if size(ph,2)==3
        ph(:,1,:,:,:)=mean(ph(:,[1,3],:,:,:),2);
    end
    
    d1b(:,1:2:end,:,:,:)=d1b(:,1:2:end,:,:,:).*repmat(exp(-1i*ph(:,1,:,:,:)),[1,ceil(size(d1b,2)/2),1,1,1]);
    d1b(:,2:2:end,:,:,:)=d1b(:,2:2:end,:,:,:).*repmat(exp(-1i*ph(:,2,:,:,:)),[1,floor(size(d1b,2)/2),1,1,1]);
     %d1bk=fft(d1b,[],1);
   d1bk=fft1c(d1b,1);