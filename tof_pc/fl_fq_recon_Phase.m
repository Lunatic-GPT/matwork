function fl_fq_recon_Phase(frecon,sig,background_mask)


im_sb=ri(frecon);

if exist('background_mask','var')
    
    a=ri(background_mask);
    
    nois=zeros(2,32);
    for i=1:32
        for j=1:2
            
            tmp=im_sb(:,:,1,j,i);
            
            nois(j,i)=std(tmp(a>0));
        end
    end
else
    nois=ones(2,32);
end

if sig>0

    im_sb=spatial_smoothing(im_sb,[sig,sig,0]);

end



%{ 
%%calculate angle first
ph_sb=angle(im_sb(:,:,:,1,:).*conj(im_sb(:,:,:,2,:)))*180/pi;

weight=mean(abs(im_sb),4);


ph_sb=mean(ph_sb.*weight,5)./mean(weight,5);
%}

%% calculate average of complex conjugate; 
%%{
res=im_sb(:,:,:,1,:).*conj(im_sb(:,:,:,2,:));

if size(im_sb,5)==32
    nois=reshape(mean(nois,1),[1,1,1,1,32]);
    sz=size(res);
    sz(5)=1;
    weight=1./repmat(nois,sz);%./mean(abs(im_sb),4);
    
    res=mean(res.*weight,5)./mean(weight,5);
end

ph_sb=angle(res)*180/pi;



prefix=strtok(frecon,'.');

sig_str=strrep(num2str(sig),'.','_');

save([prefix,'_Phase_sig',sig_str],'ph_sb');
%}

%% calculate signal average across coils first; This method is not good;
%{
res=mean(im_sb(:,:,:,1,:),5).*conj(mean(im_sb(:,:,:,2,:),5));

ph_sb=angle(res)*180/pi;

prefix=strtok(frecon,'.');

sig_str=strrep(num2str(sig),'.','_');

save([prefix,'_AveSgnl_sig',sig_str],'ph_sb');

%}



