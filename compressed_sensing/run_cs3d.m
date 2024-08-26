function im_res=run_cs3d(z,z_init,z_true,ind_sparse,mask,XFM,TVWeight,xfmWeight,niter,nodisplay)

%im_res=run_cs(z,z_init,z_true,ind_sparse,mask,XFM,TVWeight,xfmWeight,niter,nodisplay)

if ~exist('TVWeight','var')
 TVWeight = 0.002; 	% Weight for TV penalty
end
if ~exist('xfmWeight','var')
 xfmWeight = 0.003;	% Weight for Transform L1 penalty
end

if ~exist('nodisplay','var')
    nodisplay=false;
end

if ~exist('niter','var')
    niter=5;
end

N = size(z); 	% image Size

%generate Fourier sampling operator
FT = p3DFTxp(mask);
data=z.*mask;

% scale data
FT2 = p3DFTxp(ones(size(mask)));
im_init=FT2'*z_init;

%tmp = FT2'*(tmp2);

hf=figure;

if ~isempty(z_true)
  res_true=XFM*(FT2'*z_true);
  % do iterations
  [tmp,ind_true]=sort(abs(res_true(:)),'Descend');
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
  plot(ord_true,abs(res_true(ind_sparse)),'k.','MarkerSize',6);
  set(gca,'xscale','log');
  set(gca,'yscale','log');

  xlim([1,length(ind_old)]);
  %ylim([0.001,20]);
end

%generate transform operator
%XFM = Wavelet_rect('Daubechies',4,4);	% Wavelet
% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOPxp;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = 10;
%param.TVtWeight=TVWeight;  
%param.TVt= TVOPt;
isl=round(size(im_init,3)/2);
if ~nodisplay
   figure(101);subplot(1,2,1); imshow(abs(im_init(:,:,isl)),[]);
  subplot(1,2,2); imshow(angle(im_init(:,:,isl)),[]); drawnow;
end

res = XFM*im_init;

    
for n=1:niter
	res = fnlCg_xp(res,param);
	im_res = XFM'*res;
    if n==100
        ph=angle(im_res);
        param.FT=p2DFTxp(mask,N,cos(ph)+1i*sin(ph));
        im_res=im_res.*(cos(ph)-1i*sin(ph));
        res=XFM*im_res;
    end
    
    if ~nodisplay
	  figure(101);subplot(1,2,1); imshow(abs(im_res(:,:,isl)),[]);
      subplot(1,2,2);imshow(angle(im_res(:,:,isl)),[]); drawnow;
    end
 %% output some information   
    zres=FT*im_res;
    [tmp,ind]=sort(abs(res(:)),'Descend');
   
    if ~isempty(z_true)
     norm1=sqrt(mean(abs(zres(:)-z_true(:)).^2));
     norm2=mean(abs(res(ind_sparse)-res_true(ind_sparse))./abs(res_true(ind_sparse)));
    
     fprintf('l2 norm of error wrt noiseless image = %f\n',norm1);
     fprintf('mean fractional error at sparse indices = %f\n',norm2);
    end
    
    ord_fit=zeros(1,length(ind_sparse));
    for k=1:length(ind_sparse)
        ord_fit(k)=find(ind_sparse(k)==ind);
    %    fprintf('COI%d: order (fit/true) = %d/%d, coef (fit/true)= %4.2f/%4.2f\n',k, ord_fit,ord_true(k),res(ind_sparse(k)),res_true(ind_sparse(k)));
    end
    
 %  mism=(ind(1:100000)~=ind_old(1:100000));
%   fprintf('number of mismatch = %d/%d\n',sum(mism),length(mism));
   
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
   
   h4=plot(ord_fit,abs(res(ind_sparse)),'g.','MarkerSize',6);drawnow;
end





