


N=[64,64];
    L=2;
    dwtmode('per'); 

    Xtmp=eye(N);

    W=zeros(N);
    for p=1:N(1)
        [W(:,p),L2] = wavedec(Xtmp(:,p),L,'db4');
    end

    %generate wavelet transform operator
    XFM = Wavelet(W);	
    
    
gamma=0.01;

dtemp=reshape(data_3d,[64,64,1,90]);
FT=p2DFTxp(ones(64,64,1,90));
y=FT*dtemp;

modcsresCausalDetection_xp(mask2,XFM,gamma,y)
 
 