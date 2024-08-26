function recon=ktFOCUSS_KLT_new(y, prev_recon,y_true,ind_sparse, M,factor,lambda, Minner, Mouter, pred)
%recon=KTFOCUSS_KLT_new(y, prev_recon,y_true,ind_sparse, M,factor,lambda,
%Minner, Mouter, pred)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Jong's k-t FOCUSS KLT code for cartesian trajectory %%%%%%%%%%%%%%%%%%%%%
% y: zero-padded measurements on k-t space, dimension=(ky, kx, time); 
% prev_recon:  reconstruction from previous iteration  (ky,kx, time);
% M: sampling mask, same size with y, sampled point:1, otherwise: 0;
%    
% factor: weighting matrix update power (0.5, 1), (default=0.5);
% lambda: regularization factor (default=0);
%         when noise is severe, increase lambda
% Minner: Iteration number of Conjugate Gradient (default=20);
% Mouter: FOCUSS iteration number (default=2);
% pred :  flag for subtracting temporal DC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if pred==1 && ~isempty(y_true)
    ymean=mean_sampled_kt(y,M);
    sz=size(y);
    sz(1:end-1)=1;
    y_true=y_true-repmat(ymean,sz);
end

y=circshift(squeeze(y),[size(y,1)/2,size(y,2)/2,0]);
y=permute(y,[2,1,3]);
M=circshift(squeeze(M),[size(y,1)/2,size(y,2)/2,0]);
M=permute(M,[2,1,3]);

prev_recon=squeeze(prev_recon);

%%
A = @(x,mask)  (fft(x,[],1).*mask);
AT = @(x,mask) (ifft(x.*mask,[],1));

%% dimension setting
if length(size(y))==3
    [ny, nx, nframe]=size(y);
    t_dim = length(size(y));
end
y=ifft(y,[],2); %
y=y.*M;



%%%%%%%%%%%% calculating temporal average image %%%%%%%%%%%
if pred == 1
    DC=sum(y,3);
    SM=sum(M,3);
    SM(find(SM==0))=1;
    DC=DC./SM;

    for i=1:nframe
        y(:,:,i)=y(:,:,i)-DC;
    end
    y=y.*M;
    DC=ifft(DC,[],1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% KLT derivation  (fft--> KLTH,  ifft --> KLT)
%{
mean_img = sum(prev_recon,nd);
for i=1:nframe
    prev_recon(:,:,i) = prev_recon(:,:,i) - mean_img;
end

img = reshape(prev_recon, ny*nx,nframe);
[U,a,b]= svd(img'*img);
%}
U=KLT_matrix(prev_recon);
%{ 
use it to check the dominant component
disp(size(U))
figure;plot(abs(U(:,1)));

%}
%% plot results added by xz.

XFM=KLT(prev_recon);
FT2=p2DFTxp(ones(size(y)));
if ~isempty(y_true)
    hf=figure;
    sz=size(y_true);
  y_truec=y_true;
  x_truec=XFM*(FT2'*y_truec);
  % do iterations
  [tmp,ind_true]=sort(abs(x_truec(:)),'Descend');
  ord_true=zeros(1,length(ind_sparse));
  for k=1:length(ind_sparse)
    ord_true(k)=find(ind_sparse(k)==ind_true);
  end

  ind_old=ind_true;


  subplot(2,1,1);
%{
plot(ind,'r.','MarkerSize',2);
set(gca,'xscale','log');
set(gca,'yscale','log');
   xlim([0,length(ind)]);
%}
  subplot(2,1,2);
  plot(tmp,'r');
  hold on;
  plot(ord_true,abs(x_truec(ind_sparse)),'k.','MarkerSize',6);
  set(gca,'xscale','log');
  set(gca,'yscale','log');

  xlim([1,length(ind_old)]);

end

%%  weighting factor initialization
W = ones(size(y));

%%%%%%%%%%% FOCUSS update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ktFOCUSS KLT ing...')
for outer=1:Mouter
    phi=zeros(size(W));
    x=zeros(size(W));
    x_old=zeros(size(W));
    for inner=1:Minner
        x = KLT_localh(x, U);
        x=A(x,M);
        r=y-M.*x;
        g = AT(r,M);
        g= KLT_local(g,U);

        g=g*ny;
        g=-W.*g + lambda*phi;
        
        %% conjugate gradient
        if (inner>1)
            beta = g(:)'*g(:)/(prev_g(:)'*prev_g(:));
            d = -g+beta*prev_d;
        else
            d= -g;
        end

        prev_d = d;
        prev_g = g;
        %%%%%%%%%%%%%%%%%%%%%%
        
        %%%% finding alpha
        z=W.*d;
        z = KLT_localh(z,U);
        z = A(z,M);
        alpha=(real(r(:)'*z(:)-lambda*(d(:)'*phi(:))))/(z(:)'*z(:)+lambda*d(:)'*d(:));
        phi=phi+alpha*d;
        x=W.*phi;
        
        
        frac_dx=abs(1-abs(x_old./x));
        frac_dx=frac_dx(~isnan(frac_dx));
        
        fprintf('CG: relative change = %f\n',mean(frac_dx));
        x_old=x;
        if mean(frac_dx)<1e-4
          %  break;
        end
    end
%%
    tmp=KLT_localh(x,U);
    %tmp=fft(x,[],t_dim);
    tmp=fft(tmp,[],2);
    y_fit=fft(tmp,[],1);
    sz=size(y_fit);
    y_fitc=circshift(y_fit,[sz(1)/2,sz(2)/2,0,0]);   
    y_fitc=permute(y_fitc,[2,1,3]);
    x_fitc=XFM*(FT2'*y_fitc);
  
    [tmp,ind]=sort(abs(x_fitc(:)),'Descend');
   
    if ~isempty(y_true)
     norm1=sqrt(mean(abs(y_fitc(:)-y_truec(:)).^2));
     
     norm2=mean(abs(x_fitc(ind_sparse)-x_truec(ind_sparse))./abs(x_truec(ind_sparse)));
    
     fprintf('l2 norm of error wrt noiseless image = %f\n',norm1);
     fprintf('mean fractional error at sparse indices = %f\n',norm2);
    end
   
    ord_fit=zeros(1,length(ind_sparse));
    for k=1:length(ind_sparse)
        ord_fit(k)=find(ind_sparse(k)==ind);
    %    fprintf('COI%d: order (fit/true) = %d/%d, coef (fit/true)= %4.2f/%4.2f\n',k, ord_fit,ord_true(k),res(ind_sparse(k)),res_true(ind_sparse(k)));
    end
    
   frac=zeros(1,50);
   r=logspace(3,log10(length(ind)),50);
   r=floor(r);
   if exist('ind_old','var')
    for j=1:length(r)  
     comm=intersect(ind(1:r(j)),ind_old(1:r(j)));
     frac(j)=length(comm)/r(j);
    end
   end
   ind_old=ind;
   if exist('h2','var')
       delete(h2);
       delete(h3);
        delete(h4);
   end
   if ~isempty(y_true)
   figure(hf);
   subplot(2,1,1);
   hold on;
   h2=plot(r,frac,'-b','MarkerSize',2); 
   set(gca,'xscale','log');
   
   xlim([1,length(ind)]);
   drawnow;
   subplot(2,1,2);
   hold on;
   set(gca,'xscale','log');
   set(gca,'yscale','log');

   h3=plot(tmp,'b','MarkerSize',2); 
   
   h4=plot(ord_fit,abs(x_fitc(ind_sparse)),'g.','MarkerSize',6);drawnow;
   end
   %%
    W=abs(x).^factor;    %%%%%%%% weighting matrix update
    W=W/max(W(:));       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['iteration #:', int2str(outer)])
end


recon = KLT_localh(x,U);

if pred == 1
    for i=1:nframe
        recon(:,:,i)=recon(:,:,i)+DC;
    end
end

%% added by XZ
recon=circshift(recon,[size(y,1)/2,size(y,2)/2,0]);
recon=permute(recon,[2,1,3]);
recon=recon*sqrt(size(recon,1)*size(recon,2));


function c= KLT_local(x,U)
%% KLT transform
[nframe,a]= size(U);
 c= zeros(size(x));
 [ny,nx,nt]=size(x); 
 tmp = reshape(x,ny*nx,nt);
 tmp = tmp*U;
 c = reshape(tmp, ny,nx,nt);

 
 
 
 
function c= KLT_localh(x,U)
%% KLT adjointtransform
[nframe,a]= size(U);
U = U';
 c= zeros(size(x));
 [ny,nx,nt]=size(x); 
 tmp = reshape(x,ny*nx,nt);
 tmp = tmp*U;
 c = reshape(tmp, ny,nx,nt);

 
 
 

    
