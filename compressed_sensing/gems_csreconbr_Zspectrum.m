function gems_csreconbr_Zspectrum(freq,kdata,m2,sigma,B0,dict,ktrue)
%gems_csreconbr_Zspectrum(freq,kdata,m,sigma,B0,dict)
%dict is a structure containing 
% D: is the n*m dictionary
% frf: 1*n frequencies
% By default, use theoretical dictionary
% B0 are the frequency shifts


tstart=tic;

%a=fov_shift_kspace(a3,[0,-of1],fov);


kdata=kdata.*m2;

mtrx=size(kdata);
%%
if ~exist('dict','var')
if exist('c:/users/xiaopeng/Dropbox/simulation/CS_CEST_l1/Dict_501pt_1DS.mat','var')
    dict=load('c:/users/xiaopeng/Dropbox/simulation/CS_CEST_l1/Dict_501pt_1DS.mat');
else
    dict=load('e:/Dropbox/simulation/CS_CEST_l1/Dict_501pt_1DS.mat');
end
    
end

if isa(B0,'char')
 B0=BrikLoad(B0);
end


D=zeros(length(freq),size(dict.D,2),mtrx(1),mtrx(2));

for B1iter=1:2
for i=1:mtrx(1)
    for j=1:mtrx(2)
        D(:,:,i,j)=interp1(dict.frf,dict.D,freq-B0(i,j));
    end
end

A=@(z) omp_mat(z,D,m2);
At=@(z) omp_matt(z,D,m2);

 Ac = @(z) counter( @(z) omp_mat(z,D,m2), z);
 Atc= @(z) counter( @(z)omp_matt(z,D,m2), z);
    
sz=size(m2);
%sigma=1;
n=prod(sz(1:3))*size(D,2);
Lambda=sigma*sqrt(2*log(n));
  La = my_normest( A, At, n)^2;
muf = 0.1*sigma/La; %--- can be chosen to be small 
% the factor was 0.1
muf=0.01*sigma/La;
opts = [];
opts.MaxIntIter = 6;
opts.TOlVar = 1e-6;
opts.verbose = 50;
opts.maxiter = 3000;

opts.U=speye(prod(sz(1:3))*size(D,2));
opts.Ut=speye(prod(sz(1:3))*size(D,2));
opts.stoptest = 1;  


counter();
if exist('ktrue','var')
 errFunc =  @(x) error_func(x,D,ktrue,m2);
 opts.errFcn=errFunc;
end

z3=kdata(m2>0);
tic;
    
[x_ref,niter,residuals,outputData,opts] =NESTA_UP(Ac,Atc,z3(:),Lambda,La,muf,opts);


N2 = counter();fprintf('Took %d calls in %f s\n',N2,toc);

x_ref=reshape(x_ref,[sz(1:3),size(D,2)]);
img=Dmult(D,x_ref);


B0=get_B0map(freq,abs(img));
B0(abs(B0)>400)=0;

save(sprintf('gems_csreconbr_Zspectrum%d',B1iter));
end


fprintf('gems_csreconbr_Zspectrum finished in %4.1f s\n',toc(tstart));
    
 function B1=get_B0map(frf,img)
 
     B1=zeros(size(img,1),size(img,2),size(img,3));
     d=img; 
     mask=(d(:,:,:,1)>max(d(:))*0.02)|(d(:,:,:,end)>max(d(:))*0.02);
      
     for i=1:size(img,1)
         for j=1:size(img,2)
             for k=1:size(img,3)
                 if mask(i,j,k)==0
                     continue;
                 end
              y4=squeeze(img(i,j,k,:));
              
              [tmp,ind]=min(y4);
              ind1=ind-6;
              ind2=ind+6;
              if ind1<=0
                  ind1=1;
              end
              if ind2>length(frf)
                  ind2=length(frf);
              end
              
              B1(i,j,k)=WASSR_fit_1d(frf(ind1:ind2),y4(ind1:ind2));
              disp([i,j,k]);
             end
         end
     end


        
        

function out=omp_mat(b,D,m)
        % b in real space
        % out in k space
     
        sz=size(m);
        b=reshape(b,[sz(1:3),size(D,2)]);
  
        c= Dmult(D,b);
        
        fc=ifft2c(c);
        
        
        out=fc(m>0); 
        out=out(:);
        
function out=omp_matt(b,D,m)
     % b in k space
     % out in real space
    
        sz=size(m);
        d=zeros(size(m));
        
      %  b=reshape(b,[sz(1:3),size(D,1)]);
        d(m>0)=b;
        fd=fft2c(d);
        
        
        out= Dmult2(D,fd);
        out=out(:);        
  
        
function out=Dmult(D,b)
     % contract the second dimension of D and the last dimension of b
     % D = n*m*matrx
     % b = mtrx*m
     % out = mtrx*n
    
     sz=size(b);
        N=prod(sz(1:3));
        out=zeros([N,size(D,1)]);
        b=reshape(b,[N,size(b,4)]);
         for i=1:N
                   out(i,:)=D(:,:,i)*transpose(b(i,:));
         end
         out=reshape(out,[sz(1:3),size(D,1)]);
         
 function out=Dmult2(D,b)
     % contract the first dimension of D and the first dimension of b
     % D = m*n*matrx
     % b = mtrx*m
     % out = mtrx*n
     
        sz=size(b);
        out=zeros([prod(sz(1:3)),size(D,2)]);
        b=reshape(b,[prod(sz(1:3)),size(b,4)]);
        

         for i=1:prod(sz(1:3))
                 out(i,:)=b(i,:)*D(:,:,i);
         end
         

         out=reshape(out,[sz(1:3),size(D,2)]);
         
function y=error_func(x,D,z_true,m)
           
           sz=size(z_true);
           x=reshape(x,[sz(1:3),size(D,2)]);
           img=Dmult(D,x);
            
           z=ifft2c(img);
           y=norm(z_true(m>0)-z(m>0))/sqrt(length(z_true(m>0)));
           
           
           
           