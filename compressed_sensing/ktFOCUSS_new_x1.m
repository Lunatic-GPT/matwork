function recon=ktFOCUSS_new_x1(y, y_true,ind_sparse,num_low_phase, M,factor,lambda, Minner, Mouter, pred)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Hong's k-t FOCUSS code for cartesian trajectory %%%%%%%%%%%%%%%%%%%%%
% recon=ktFOCUSS_new(y, y_true,ind_sparse,num_low_phase, M,factor,lambda,
% Minner, Mouter, pred)
% y: zero-padded measurements on k-t space, dimension=(ky, kx, time);
%    Assumed low frequency is positioned at vertex of k-space;
%    1<=ky<=4, end-3<=ky<=end : fully sampled region 
%    ky : phase encoding 
%    kx : frequency encoding
%    time : time frames
% 
% M: sampling mask, same size with y, sampled point:1, otherwise: 0;
%    SAMPLING PATTERN should be RANDOM!!!
%    1<=ky<=4, end-3<=ky<=end : M(ky, kx, time)=1 
%    
% factor: weighting matrix update power (0.5, 1), (default=0.5);
% lambda: regularization factor (default=0);
%         when noise is severe, increase lambda
% Minner: Iteration number of Conjugate Gradient (default=20);
% Mouter: FOCUSS iteration number (default=2);
% pred  : use prediction with temporal DC
% num_low_phase:  low frequency range for W initialization
% y_true: same format as y. [] to skip
% ind_sparse: indices defined for centered images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_fit=false;


y=circshift(squeeze(y),[1,size(y,2)/2,0]);
y=permute(y,[2,1,3]);
M=circshift(squeeze(M),[1,size(M,2)/2,0]);
M=permute(M,[2,1,3]);


%% dimension setup
if length(size(y))==3
    [ny, nx, nframe]=size(y);
    t_dim = length(size(y));
end
y=ifft(y,[],2); %
y=y.*M;


A = @(x,mask)  (fft(x,[],1).*mask);
AT = @(x,mask) (ifft(x.*mask,[],1));

%%%%%%%%%%%% calculating temporal average image %%%%%%%%%%%
if pred == 1
    DC=sum(y,3);
    SM=sum(M,3);
    SM(find(SM==0))=1;
    DC=DC./SM;

    for i=1:nframe
        y(:,:,i)=y(:,:,i)-DC;
        if ~isempty(y_true);
     %     y_true(:,:,i)=y_true(:,:,i)-fft(DC,[],2);
       
        end
    end
    
    y=y.*M;
    DC=ifft(DC,[],1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot results added by xz.
if show_fit
hf=figure;
XFM=FTt;
FT2=p2DFTxp(ones(size(y)));
if ~isempty(y_true)
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
end

%%%%%%%%%%%% Initializing weighting matrix %%%%%%%%%%%%%%%%
W=zeros(size(y));
offset = floor(num_low_phase/2);
W(1:offset+1,:,:)=y(1:offset+1,:,:);
W(end-offset:end,:,:)=y(end-offset:end,:,:);
W=ifft(W,[],1);
W=ifft(W,[],3);
W=abs(W).^factor;
% if pred == 1
%     W(:,:,1)=0;
% end
W=W/max(W(:));
%W=ones(size(y));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% W= abs(ifft(AT(low_resol_y,M),[],t_dim)).^factor;
% W=W/max(W(:));


%%%%%%%%%%% FOCUSS update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ktFOCUSS ing...')

for outer=1:Mouter
    phi=zeros(size(W));
    
    x=zeros(size(W));
    x_old=zeros(size(W));
    
    for inner=1:Minner
        x=fft(x,[],t_dim);
        x=A(x,M);
        r=y-M.*x;
        g = AT(r,M);
        g=ifft(g,[],t_dim);
        
        g=g*nframe*ny;   %% this is very important !!! by Jong.
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
        z=fft(z,[],t_dim);
        z = A(z,M);
        alpha=(real(r(:)'*z(:)-lambda*(d(:)'*phi(:))))/(z(:)'*z(:)+lambda*d(:)'*d(:));
        phi=phi+alpha*d;
        x=W.*phi;
        
        frac_dx=abs(1-abs(x_old./x));
        frac_dx=frac_dx(~isnan(frac_dx));
        
        fprintf('%d-CG: relative change = %f\n',inner,mean(frac_dx));
        
        x_old=x;
        if mean(frac_dx)<1e-5
        %    break;
        end
    end
 %% plot    added by xz
 if show_fit
    tmp=fft(x,[],t_dim);
    tmp=fft(tmp,[],2);
    y_fit=fft(tmp,[],1);
    
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
 
 %% end plot  
    W=abs(x).^factor;    %%%%%%%% weighting matrix update
    W=W/max(W(:));       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['iteration #:', int2str(outer)])
end



recon=fft(x,[],t_dim);
if pred == 1
    for i=1:nframe
        recon(:,:,i)=recon(:,:,i)+DC;
    end
end

%% added by XZ
recon=circshift(recon,[size(y,1)/2,size(y,2)/2,0]);
recon=permute(recon,[2,1,3]);
recon=recon*sqrt(size(recon,1)*size(recon,2));



        