function gems_csreconbr_Zspectrum_dct(kdata,m2,sigma,ktrue,fname)
%gems_csreconbr_Zspectrum(kdata,m2,sigma,ktrue)
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


nt=size(kdata,4);

D=@(x) idct_trunc(x,1);
Dt=  @(x) dct_trunc(x,1);

%D=@(x) iWavelet1D(x);
%Dt=  @(x) Wavelet1D(x);

A=@(z) omp_mat(z,D,m2);
At=@(z) omp_matt(z,Dt,m2);

 Ac = @(z) counter( @(z) omp_mat(z,D,m2), z);
 Atc= @(z) counter( @(z)omp_matt(z,Dt,m2), z);
    
sz=size(m2);
%sigma=1;
n=prod(sz(1:3))*nt;
Lambda=sigma*sqrt(2*log(n));
  La = my_normest( A, At, n)^2;
muf = 0.1*sigma/La; %--- can be chosen to be small 
% the factor was 0.1
muf=0.01*sigma/La;
opts = [];
opts.MaxIntIter = 6;
opts.TOlVar = 1e-6;
opts.verbose = 10;
opts.maxiter = 3000;

opts.U=speye(prod(sz(1:3))*nt);
opts.Ut=speye(prod(sz(1:3))*nt);
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

x_ref=reshape(x_ref,[sz(1:3),nt]);
img=Dmult(D,x_ref);

write_afni(abs(img),fname);


fprintf('gems_csreconbr_Zspectrum finished in %4.1f s\n',toc(tstart));
    

function out=omp_mat(b,D,m)
        % b in real space
        % out in k space
     
        sz=size(m);
        b=reshape(b,[sz(1:3),size(m,4)]);
  
        c= Dmult(D,b);
        
        fc=ifft2c(c);
        
        
        out=fc(m>0); 
        out=out(:);
        
function out=omp_matt(b,D,m)
     % b in k space
     % out in real space
    
        d=zeros(size(m));
        
      %  b=reshape(b,[sz(1:3),size(D,1)]);
        d(m>0)=b;
        fd=fft2c(d);
        
        
        out= Dmult(D,fd);
        out=out(:);        
  
        
function out=Dmult(D,b)
     % contract the second dimension of D and the last dimension of b
     % D = n*m*matrx
     % b = mtrx*m
     % out = mtrx*n
    
     sz=size(b);
        N=prod(sz(1:3));
        out=zeros([N,size(b,4)]);
        b=reshape(b,[N,size(b,4)]);
         for i=1:N
                   out(i,:)=D(b(i,:));
         end
         
         out=reshape(out,sz);
         
function y=error_func(x,D,z_true,m)
           
           sz=size(z_true);
           x=reshape(x,[sz(1:3),size(z_true,4)]);
           img=Dmult(D,x);
            
           z=ifft2c(img);
           y=norm(z_true(m>0)-z(m>0))/sqrt(length(z_true(m>0)));
           
           
           
           