function res=rCBF_pCASL_Artery(dname,nslice,mosaic)
%% mosaic: [rows,cols]

if ~exist('mosaic','var')
    mosaic=[1,1];
end

if ~exist('nslice','var')
    nslice=1;
end


    d=ri(dname,nslice);
    sl=d(:,:,:,1:2:end);
    
    sc=d(:,:,:,2:2:end);
    
    nd=ndims(sc);
    res=mean((sc-sl)./sc,nd);

    res=permute(res,[2,1,3]);

    sz=size(res);    
   cols=mosaic(2);
   rows=mosaic(1);
   
    res=reshape(res,[sz(1)/cols,cols,sz(2)/rows,rows,sz(3:end)]);
    
    res=permute(res,[1,3,2,4,5]);
    res=reshape(res,[sz(1)/cols,sz(2)/rows,length(res(:))/sz(1)/sz(2)*cols*rows]);
    save([dname,'_CBF'],'res');
    
    
    
    