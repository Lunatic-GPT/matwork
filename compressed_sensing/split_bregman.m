%           mrics.m  by Tom Goldstein (TomGoldstein1@gmail.com)
%     This file contains methods for performing compressed sensing
%  recontructions of images from k-space data using the Split Bregman 
%  method.  
%     To use the method, simply add this "m" file to your current directory, 
%  and then call the following method:
%   
%              u = mrics(R,F, mu, lambda, lambdat,gamma, nInner, nOuter,xfm,c1);
%   
%  The inputs to this function are described below:
%  
%  R - This is a matrix which determines which elements of K-Space
%           are known.  Take R(m,n)=0 if Fourier mode (m,n) is unknown,
%           and R(m,n)=1 if the corresponding mode is known.
%  F - This is the K-Space (Fourier) data.  In other words, F is the 
%           Fourier transform of the image you wish to recover.  If a
%           Fourier mode is known, then it should have a non-zero value.
%           If a Fourier mode is unknown, then simply set the corresponding
%           entry in this matrix to zero.  If you have set the values
%           in this matrix properly, then you should have (R.*F==F).
%  mu- The parameter on the fidelity term in the Split Bregman method.  The
%           best choice will depend on how the data is scaled.
% lambda - The coefficient of the constraint term in the Split Bregman
%      model. For most problems, I suggest using lambda=mu.
% gamma - This is a regularization parameter.  I suggest that you take
%      gamma = mu/100.
% nInner - This determines how many "inner" loops the Split Bregman method
%      performs (i.e. loop to enforce the constraint term).  I suggest
%      using nInner = 30 to be safe.  This will usually guarantee good
%      convergence, but will make things a bit slow.  You may find that you
%      can get away with nInner = 5-10
% nOuter - The number of outer (fidelity term) Bregman Iterations.  This
%      parameter depends on how noisy your data is, but I find that
%      nOuter=5 is usually about right.
%xfm - the wavelet transform
%c1 - weighting factor for xfm.

function u = split_bregman(R,f, mu, lambda, gamma, nInner, nBreg,xfm,c1,phase)
    
if ~exist('phase','var')
    phase=1;
end
    f=f.*R;
    

    sz=size(f);
    rows=sz(1);
    cols=sz(2);
    f0 = f;
    
    
in=sum(f,4);
nm=sum(R,4);
nm(nm==0)=1;
in=in./nm;

    fftxfm=p2DFTxp(ones(size(f)),size(f),phase);
    
   u=fftxfm'*repmat(in,[1,1,1,size(f,4)]);

  %  u = zeros(sz);
    x = zeros(sz);
    y = zeros(sz);
    z = zeros(sz);
    t = zeros(sz);
    
    w=zeros(sz);
    bx = zeros(sz);
    by = zeros(sz);
    bz = zeros(sz);
    bt = zeros(sz);
    bw=zeros(sz);

     % Build Kernels
        %ph=angle(im_res);
    %fftxfm=fft2o; 

    murf = fftxfm'*(mu*R.*f);
    
   % uker = zeros(rows,cols);
   % uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
    % the above matrix is same as: dx'*dx*ifft2(ones(N,N))+ifft2(ones(N,N))*dy*dy'
    if lambda(1)>0
     ukerxy=Dxt(Dx(fftxfm'*ones(size(u))))+Dyt(Dy(fftxfm'*ones(size(u))));
     uker=mu*R+lambda(1)*(fftxfm*ukerxy)+gamma+c1;
    else
     uker=mu*R+gamma+c1;
    end    
        
    if length(lambda)>=2&&lambda(2)>0
     ukerz=Dzt(Dz(fftxfm'*ones(size(u))))*(lambda(2)/lambda(1))^2;
     ukerz=fftxfm*ukerz;
     uker=uker+ukerz*lambda(1);
    end
    if length(lambda)==3&&lambda(3)>0
     ukert=Dtt(Dt(fftxfm'*ones(size(u))));
     ukert=fftxfm*ukert;
     uker=uker+lambda(3)*ukert;
    end
    
   %%initial 
   tv=0;
        if lambda(1)>0
            dx = Dx(u);
            dy =Dy(u);
            tv=sqrt(abs(dx).^2+abs(dy).^2);
            tv=sum(tv(:));
        end
        tvz=0;
        if length(lambda)>=2&&lambda(2)>0
            dz=Dz(u)*lambda(2)/lambda(1);
            tvz=abs(dz);
          tvz=sum(tvz(:));
        end
        
        if c1>0
            wvl=xfm*u;
            wvl=sum(abs(wvl(:)));
        else
            wvl=0;
        end
        
        resid=(f0-R.*(fftxfm*u)).^2;
        resid=sum(abs(resid(:)));
        rms=sqrt(resid/length(find(R>0)));
        fprintf('0 s: total = %f, tv = %f, tvz = %f, data obj = %f, rms = %f, xfm= %f, \n',resid+(wvl+tv)/mu,tv,tvz,resid,rms,wvl);    
        
    %  Do the reconstruction
    for outer = 1:nBreg;
       
        for inner = 1:nInner;
             % update u
        %     disp([outer,inner]);
              tic;
             rhs = murf++gamma*u;
             
             if lambda(1)>0
                 rhs=rhs+lambda(1)*Dxt(x-bx)+lambda(1)*Dyt(y-by);
             end
             if  length(lambda)==3&&lambda(3)>0 
              rhs = rhs+lambda(3)*Dtt(t-bt);
              dt=Dt(u);
             end
             if c1>0
              rhs = rhs+c1*(xfm'*(w-bw));           
             end
             if length(lambda)>=2&&lambda(2)>0
                 rhs=rhs+lambda(2)*Dzt(z-bz);
             end
             
             if (inner>1||outer>1)
               u = fftxfm'*(fftxfm*rhs./uker);
             end
            if lambda(1)>0
            dx = Dx(u);
            dy =Dy(u);
            end
            
            if length(lambda)>=2&&lambda(2)>0
                 dz=Dz(u)*lambda(2)/lambda(1);
             end

           %fprintf('(%d, %d): rms error = %f\n',outer, inner,rms);
      %     dxrms2=lambda/2*sum(abs(x(:)-dx(:)-bx(:)).^2)/length(R(:));
      %     dyrms2=lambda/2*sum(abs(y(:)-dy(:)-by(:)).^2)/length(R(:));
      %     dxyrms=sum(sqrt(dx(:)'*dx(:)+dy(:)'*dy(:)))/length(R(:));
         % fprintf('rms = %6.5f, dxyrms = %6.5f, dxrms2 = %6.5f, dyrms2 =
         % %6.5f\n',rms,dxyrms,dxrms2,dyrms2)
            % update x and y
            
            if length(lambda)==3&&lambda(3)>0
              t=shrink(dt+bt,1/lambda(3));
              bt=bt+dt-t;
            end
            if length(lambda)>=2&&lambda(2)>0
              [x,y,z]=shrink3(dx+bx,dy+by,dz+bz,1/lambda(1));
             % z=shrink(dz+bz,1/lambda(2));
              
            % [x,y] = shrink2( dx+bx, dy+by,1/lambda(1));   
              bz=bz+dz-z;
            elseif lambda(1)>0
             [x,y] = shrink2( dx+bx, dy+by,1/lambda(1)); 
             bx = bx+dx-x;
             by = by+dy-y;
            end
            
          if c1>0
            w=shrink(xfm*u+bw,1/c1); 
            bw = bw+xfm*u-w;
          end
            
        %    fprintf('data = %f; tvx = %f; tvy = %f; wvl = %f; tvxy = %f; total = %f\n',tmp1,tmp2,tmp3,tmp4,tmp5,(tmp1+tmp2+tmp3+tmp4+tmp5));

            
        %    fprintf('bx = %f; by = %f\n',sum(abs(bx(:))),sum(abs(by(:))));
       
        if lambda(1)>0
            tv=sqrt(abs(dx).^2+abs(dy).^2);
            tv=sum(tv(:));
        else
            tv=0;
        end
        
        if length(lambda)>=2&&lambda(2)>0
          tvz=abs(dz);
          tvz=sum(tvz(:));
        else    
          tvz=0; 
        end
        
        if c1>0
            wvl=xfm*u;
            wvl=sum(abs(wvl(:)));
        else
            wvl=0;
        end
        
        resid=(f0-R.*(fftxfm*u)).^2;
        resid=sum(abs(resid(:)));
        rms=resid/length(find(R>0));
        fprintf('%4.1f s: total = %f, tv = %f, tvz = %f, data obj = %f, rms = %f, xfm= %f\n',toc,resid+(wvl+tv)/mu,tv,tvz,resid,rms,wvl);      
        figure(103); subplot(1,2,1);imshow(abs(u(:,:,1,1)),[]);
        subplot(1,2,2);imshow(angle(u(:,:,1,1)),[]);
        
        drawnow;
        end
        fprintf('-----------------------------------------------------------------------------------\n');
        f = f+f0-R.*(fftxfm*u);
        murf = fftxfm'*(mu*R.*f);
       % disp(toc);
     
    end
    


function d = Dt(u)

if ndims(u)==3
    d=u-u(:,:,[end,1:end-1]);
elseif ndims(u)==4
    d=u-u(:,:,:,[end,1:end-1]);
end


function d = Dtt(u)

if ndims(u)==3
    d=u-u(:,:,[2:end,1]);
elseif ndims(u)==4
    d=u-u(:,:,:,[2:end,1]);
end


function d = Dx(u)
cols = size(u,2);
d = zeros(size(u));
d(:,2:cols,:,:) = u(:,2:cols,:,:)-u(:,1:cols-1,:,:);
d(:,1,:,:) = u(:,1,:,:)-u(:,cols,:,:);
return

function d = Dxt(u)
cols=size(u,2);
d = zeros(size(u));
d(:,1:cols-1,:,:) = u(:,1:cols-1,:,:)-u(:,2:cols,:,:);
d(:,cols,:,:) = u(:,cols,:,:)-u(:,1,:,:);
return


function d = Dy(u)
rows = size(u,1); 
d = zeros(size(u));
d(2:rows,:,:,:) = u(2:rows,:,:,:)-u(1:rows-1,:,:,:);
d(1,:,:,:) = u(1,:,:,:)-u(rows,:,:,:);
return

function d = Dyt(u)
rows = size(u,1);
d = zeros(size(u));
d(1:rows-1,:,:,:) = u(1:rows-1,:,:,:)-u(2:rows,:,:,:);
d(rows,:,:,:) = u(rows,:,:,:)-u(1,:,:,:);
return


function d = Dz(u)
cols = size(u,3);
d = zeros(size(u));
d(:,:,2:cols,:) = u(:,:,2:cols,:)-u(:,:,1:cols-1,:);
d(:,:,1,:) = u(:,:,1,:)-u(:,:,cols,:);
return

function d = Dzt(u)
cols=size(u,3);
d = zeros(size(u));
d(:,:,1:cols-1,:) = u(:,:,1:cols-1,:)-u(:,:,2:cols,:);
d(:,:,cols,:) = u(:,:,cols,:)-u(:,:,1,:);
return

function [xs,ys] = shrink2(x,y,lambda)

s = sqrt(x.*conj(x)+y.*conj(y));
ss = s-lambda;
ss = ss.*(ss>0);

s = s+(s<lambda);
ss = ss./s;

xs = ss.*x;
ys = ss.*y;

return;

function xs=shrink(x,lambda)

s=abs(x);
ss = s-lambda;
ss = ss.*(ss>0);

s = s+(s<lambda);
ss = ss./s;

xs = ss.*x;

function [xs,ys,zs] = shrink3(x,y,z,lambda)

s = sqrt(x.*conj(x)+y.*conj(y)+z.*conj(z));
ss = s-lambda;
ss = ss.*(ss>0);

s = s+(s<lambda);
ss = ss./s;

xs = ss.*x;
ys = ss.*y;
zs=ss.*z;

return;
