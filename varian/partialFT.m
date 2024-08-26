function raw=partialFT(raw,negphases)
% raw=partialFT(raw,negphases)
        % adopted from Cecil's epirecon code.
        % partial FT along the second dimension in raw
        % the results are already FT'ed along the first and second dim
np=size(raw,1)*2;
nnv=size(raw,2);
siz=size(raw);
arraydim=size(raw,4);
ns=size(raw,3);
raw=fftshift(raw,2);
        % 1 Create a truncated and zero padded k-space data(AAzf) from the shiftecho data(AAraw)
   raw = ifft(ifftshift(raw,1),np/2,1)*np/2; % swap in RO direction before fft

          AAraw = fftshift(fftshift(fft(raw,[],1),1),2)/np*2;
        %AAraw = fftshift(fft(raw,[],1),1)/np*2;
        %AAraw = fftshift(raw,2)/np*2; 
        lpe=1;
ppe=-lpe/2;        
           % fprintf('shifting ppe=%gmm...',ppe*10)
    peShift = ones(np/2,1)*exp(-1i*2*pi*ppe*(0:nnv-1)/lpe);
    for m=1:arraydim
        for n=1:ns
            raw(:,:,n,m) = raw(:,:,n,m).*peShift;
        end
    end

        
        siz(2) = (nnv-negphases)*2;
        AAzf = zeros(siz);
        cmplx_new = AAzf;
        cmplx_old = AAzf;
        AAzf(:,nnv-negphases*2+1:nnv,:,:) = AAraw(:,1:negphases*2,:,:);
        
        % 2 Generate a low freq. phase map(lowph)
        AAlowph = angle(ifft2(ifftshift(ifftshift(AAzf,1),2),np/2,siz(2)));
        
        % 3 Use the aymmetric plus zero-filling data as our inital magnitude data
        AAzf(:,nnv+1:end,:,:)=AAraw(:,negphases*2+1:end,:,:);
        
        % 4 Generate an initial complex data(cmplx)
        AAcmplx = abs(ifft2(fftshift(fftshift(AAzf,1),2),np/2,siz(2))).*exp(AAlowph*1i);
        
        for tt=1:arraydim
            for slc=1:ns
                disp([tt,slc]);
                threshold=sqrt(var(reshape(AAzf(:,:,slc,tt),size(AAzf,1)*size(AAzf,2),1)))*3;
                for itr=1:50
                    % 5 Create an intermediate k-space data
                    AAzf(:,:,slc,tt) = fft2(AAcmplx(:,:,slc,tt),np/2,siz(2));
                    
                    % 6 Subsitute part of the intermediate k-space data into shiftecho data(raw_new)
                    AAzf(:,:,slc,tt)=ifftshift(ifftshift(AAzf(:,:,slc,tt),1),2);
                    AAzf(:,nnv-negphases*2+1:end,slc,tt)=AAraw(:,:,slc,tt);
                    AAzf(:,:,slc,tt)=fftshift(fftshift(AAzf(:,:,slc,tt),1),2);
                    
                    % 7 fft to create new complex data
                    cmplx_new(:,:,slc,tt) = ifft2(AAzf(:,:,slc,tt),np/2,siz(2));
                    
                    if sum(sum(abs(cmplx_new(:,:,slc,tt)-cmplx_old(:,:,slc,tt))))<threshold
                        AAzf(:,:,slc,tt)=ifftshift(ifftshift(AAzf(:,:,slc,tt),1),2);
                        break
                    end
                    cmplx_old(:,:,slc,tt)=cmplx_new(:,:,slc,tt);
                    
                    % 10 Generate an intermediate complex data(cmplx) and iterate
                    AAcmplx(:,:,slc,tt) = abs(cmplx_new(:,:,slc,tt)).*exp(AAlowph(:,:,slc,tt)*1i);
                end
            end
        end
        
        % 9 Smoothen via uu-point Hanning filter and output the result to AAcmplx
        uu = 3;
        for pp = 1:uu
            AAzf(:,nnv-negphases*2+pp,:,:) = (AAraw(:,pp,:,:)*(1+cos(pi+pi*(pp-1)/uu)) + AAzf(:,pp,:,:)*(1+cos(pi*(pp-1)/uu)))/2;
        end
        
        hamm2d=(25-21*cos(2*pi*(0:siz(1)-1)'/siz(1)))/46*(25-21*cos(2*pi*(0:siz(2)-1)/siz(2)))/46;
       if length(siz)==2 
           siz(3)=1;
           siz(4)=1;
       end
       if length(siz)==3 
           siz(4)=1;
       end
        for p=1:siz(3)
            for q=1:siz(4)
               % AAzf(:,:,p,q)=AAzf(:,:,p,q).*hamm2d;
            end
        end
        
     %   raw = fftshift(fftshift(ifft2(fftshift(fftshift(AAzf,1),2),np/2,siz(2)),1),2)*siz(2)*np/2;
        %raw = fftshift(ifft2(fftshift(fftshift(AAzf,1),2),np/2,siz(2)),1)*siz(2)*np/2;
        
      raw = ifft2(fftshift(fftshift(AAzf,1),2),np/2,siz(2))*siz(2)*np/2;
      
%raw=raw/sqrt(size(raw,2));
%img = abs(raw(:,:,:,:));
%ph=angle(raw(:,:,:,:));


