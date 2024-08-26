function csrecon_Hongtu()
%gems_csrecon(a[,fB1,sigma])

sigma=1;


dict=load('Dict_501pt_1DS.mat');
    


waterf=load('WaterFreq.mat');
waterf=waterf.B0;    %water resonance frequences.

data=load('14_RawData.mat');
frf2=data.frequency; %offset frequencies

m2=data.m2;  %kspace undersampling mask: 1 sampled; 0 not sampled.
z2=data.kdata; % undersampled k-space data.

mtrx=[64;64];
D=zeros(length(frf2),size(dict.D,2),mtrx(1),mtrx(2));


for B1iter=1:2
for i=1:mtrx(1)
    for j=1:mtrx(2)
        D(:,:,i,j)=interp1(dict.frf,dict.D,frf2-waterf(i,j));  % theorectical dictionary
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


z3=z2(m2>0);

x_ref =NESTA_UP(Ac,Atc,z3(:),Lambda,La,muf,opts);


N2 = counter();
x_ref=reshape(x_ref,[sz(1:3),size(D,2)]);
img=Dmult(D,x_ref);
waterf=get_B1map(frf2,abs(img));
waterf(abs(waterf)>500)=0;


%write_afni(abs(img),[a,'_recon']);

end

save('recon_results');


   function B1=get_B1map(frf,img)
 
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


              
        
function d2=phaseCorr(d2,nav)

  ph=angle(nav);
  
  ph=repmat(ph,[size(d2,1),1,1,1]);
  
 d2=d2.*exp(-1i*ph);

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
           
           
           
           