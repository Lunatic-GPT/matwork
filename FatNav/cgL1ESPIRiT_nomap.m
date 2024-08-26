function [res] = cgL1ESPIRiT_nomap(kData, x0, FT,  nIterCG, XOP, lambda, alpha,nIterSplit)
%
%[res] = cgESPIRiT(kData, x0, FT, MapOp, nIterCG, [ XOP, lambda, alpha,nIterSplit)
%
% Implementation of image-domain L1-Wavelet regularized ESPIRiT reconstruction from arbitrary
% k-space. The splitting is based on F. Huang MRM 2010;64:1078?1088
% Solves the problem: || Em - y ||^2 + \lambda ||\Psi x||_1 + \alpha||x-m||^2 
% by splitting into two subproblems:
%   I) ||Em - y || ^2 + \alpha ||x-m||^2
%   II) \alpha||x-m||^2 + \lambda ||\Psi x||_1
%
% large \alpha has slower splitting iterations and faster cg ones. 
% 
% 
% 
% Inputs: 
%       kData   - k-space data matrix it is 3D corresponding to [readout,interleaves,coils] 
%       x0      - Initial estimate of the coil images
%       FT - fft/nufft operator (see @NUFFT @p2DFT classes)
%       nIterCG   - number of LSQR iterations
%       MapOp     - ESPIRiT Operator (See @ESPIRiT)
%       XOP - transform thresholding operator for l1 based soft
%             thresholding
%       lambda   - L1-Wavelet penalty
%       alpha    - splitting parameter (0.5 default)
%       nIterSplit - number of splitting iterations
%       
% Outputs:
%       res - reconstructed ESPIRiT images
%       
%
% 
%
% (c) Michael Lustig 2006, modified 2010, 2013

if nargin < 6
    XOP = 1;
    lambda = 0;
    alpha = 0;
    nIterSplit = 1;
end

%     
% N = prod(size(x0));
% M = prod(size(kData));
imSize = size(x0);
% m=abs(kData)>0;
% make dyadic size if Wavelets are used. 
if strcmp(class(XOP),'Wavelet') == 1

    if length(imSize)>2
        imSize_dyd = [max(2.^ceil(log2(imSize(1:2)))), max(2.^ceil(log2(imSize(1:2)))),imSize(3)];
    else
        imSize_dyd = [max(2.^ceil(log2(imSize(1:2)))), max(2.^ceil(log2(imSize(1:2))))];
    end
else
    imSize_dyd = imSize;
end

    
dataSize = [size(kData)];

res = double(x0(:));

for n=1:nIterSplit
    tic;
    if alpha>0
        b = double([kData(:); sqrt(alpha)*res(:)]);
    else
        b = double(kData(:));
    end
  %  [res,FLAG,RELRES,ITER,RESVEC,LSVEC] = lsqr(@(x,tflag)afun(x,FT,dataSize, imSize, alpha, tflag), b, [], nIterCG,speye(N,N),speye(N,N), res(:));
      [res,FLAG,RELRES,ITER,RESVEC,LSVEC] = lsqr(@(x,tflag)afun(x,FT,dataSize, imSize, alpha, tflag), b, [], nIterCG,[],[], res(:));
    res = reshape(res,imSize);

        
        
    if lambda > 0
        tmp = zpad(res,imSize_dyd);
        tmp = XOP*tmp;
        tmp = SoftThresh(tmp,lambda/sqrt(alpha));
        res = XOP'*tmp;
        res = reshape(res,imSize_dyd);
        res = crop(res,imSize);
    end
    
    obj1 = (FT* ( res) - kData); 
    
 %   figure(100), imshow3(abs(res),[]), drawnow;
    fprintf('Iteration: %d, remaining time = %d s; relative residual = %f, consistency (per pixel): %f (%f);\n ',n,round(toc*(nIterSplit-n)),RELRES,norm(obj1(:)),norm(obj1(:))/sqrt(sum(kData(:)>0)) );
                                                                
end



function [y, tflag] = afun(x,FT, dataSize, imSize,  alpha, tflag)

    
    if strcmp(tflag,'transp')
        
        
        y = reshape(x(1:prod(dataSize)),dataSize);
        xtmp = x(prod(dataSize)+1:end);
        
        x = FT'.*double(y);
      
        if alpha>0
        y = x(:)+ sqrt(alpha)*xtmp(:);
        else
            y=x(:);
        end
    else
        
        x = reshape(x,imSize);

        y = FT.*double(x);
        if alpha>0
            y = [y(:); sqrt(alpha) * x(:)];
        else
            y=y(:);
        end
    end
