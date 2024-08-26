function res = SAKE_nufft(DATA, k, imSize,kSize, wnRank,nIter)
%
% res = SAKE(DATA, kSize, wnRank,nIter)
%
% function performs an autocalibration reconstrution without calibration
% lines based on the SAKE method (P. Shin et. al, "Calibrationless Parallel 
% Imaging Reconstruction Based on Structured Low-Rank Matrix Completion"
% 2013, submitted to MRM. 
%
% It is recommended that this method be used to complete a calibration
% region and then use ESPIRiT to generate ESPIRiT maps. 
%
% INPUTS:
%       DATA [n X nc]  - 2D multi-coil k-space data with missing
%                               entries, preferably random! missing entries
%                               should be EXACTLY = 0!
%       k   [n x 2]
%       imSize [1 x 2]   size of the calibration region
%       kSize [ 1 x 2]       - sliding window size
%       wnRank [ 1 x 1]       - the window-normalized rank to enforce 
%                               (# of coils = full rank, 2.0 would be typical 
%                               for 8 coils with window size [5 5]
%       nIter                - number of iTerations.
%       show                 - show intermediate reconstruction show=0 will
%                               skip. Show = 100 will plot in figure 100
%
%
% Outputs:
%       res -                - 2D multi coil k-space data in a Cartesian grid and where the missing
%                               data was filled. 


%NUFFT(k,w,shift,imSize)


nufft=NUFFT(k,1,[0,0],imSize);

res = fft_nuifft(DATA,nufft);
sz=size(res);


for n=1:nIter
    
    % reorder data to get Hankel structure. 
    tmp = im2row(res,kSize); [tsx,tsy,tsz] = size(tmp);
    A = reshape(tmp,tsx,tsy*tsz);
    
    % SVD thresholding
    [U,S,V] = svd(A);
    keep = 1:floor(wnRank*prod(kSize));
    A = U(:,keep)*S(keep,keep)*V(:,keep)';
    
    % Enforce Hankel structure
    A = reshape(A,tsx,tsy,tsz);
    tmp = row2im(A,sz,kSize);
    
    % enforce data consistency
    
    resid=(nufft_ifft(tmp,nufft)-DATA);
    res = tmp - fft_nuifft(resid,nufft);
    
end

function res=fft_nuifft(kdata,nufft)

 for i=1:size(kdata,2)
   im=nufft'*kdata(:,i);
   res(:,:,i)=fft1c(fft1c(im,1),2);
 end

function res=nufft_ifft(kdata_cartesian,nufft)

   im=ifft1c(ifft1c(kdata_cartesian,1),2);
   for i=1:size(im,3)
    res(:,i)=nufft*im(:,:,i);
   end




