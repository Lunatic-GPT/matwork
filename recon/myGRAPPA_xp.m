function [Rk,coef] = myGRAPPA_xp(Sk,Slines,acs,kernelY,kernelX,coef,do_recon,plot_fit)

%myGRAPPA
% 
%      [Rk,coef] = myGRAPPA_xp(Sk,Slines,acs,kernelY,kernelX,coef[,do_recon,plot_fit])
%  
% do_recon:  if true (default), the missing lines will be calculated. 
% Otherwise only calculate coef.  When false, coef should be 0.
% INPUTS:
%	Sk: [nRO,nPE,coils] Original size of 2D slice
%       
%       Slines:  Position of sampled lines along PE; 1 based
%       ACS: acs lines; 1 based;
%       kernelY: kernel along PE
%       kernelX: kernel along RO
%       coefs: GRAPPA coefficientes (coef=0, then estimated)
%
% OUTPUT
%	Rk: Reconstructed k-space
%	Rx: reconstructed (complex) x-space
%	coef: GRAPPA reconstruction coefficients
%	
% GRAPPA reconsrtuction of subsampled k_space
%
%
% Implementation of  
% 
%     M. A. Griswold, P. M. Jakob, et. al..
%     Generalized autocalibrating partially parallel acquisitions
%     Mag. Reson. Med.   47(6):1202-1210.  Jun 2002
%
% uses 3x2 kernel 
% REQUIRES bicg.m (from sparfun)
%
% EXAMPLE
%       [k_rG2,Ix_G,G_coefs]=myGRAPPA([Mx My coils],k_sG,Samp,ACS,0);
%
% Based on "recongrappa.m" by Scott Hoge  (shoge at ieee dot org)
%
% PARALLEL MRI TOOLBOX
%
% Santiago Aja-Fernandez, LPI
% www.lpi.tel.uva.es/~santi
% Valladolid 23/06/2011
%


if ~exist('plot_fit','var')
    plot_fit=false;
end

if ~exist('do_recon','var')
    do_recon=true;
end

if ~do_recon  && coef~=0
    error('do_recon should be true when coef is not zero');
end

Sk=permute(Sk,[2,1,3]);
lx=length(kernelX);  %along RO
ly=length(kernelY);  %along PE

  Mx=size(Sk,1);  % nPE
  My=size(Sk,2);  % nRO
  Ncoils=size(Sk,3);

%IF acs adjustent needed----
%if rem(rate,2)
%acs=acs(1:end-1);
%else
%acs=acs(2:end-1);
%end
%----------------------------


ImK = zeros( size(Sk,1), 1 );         
ImK(Slines) = 1;  %Sampled lines set to 1

acs2=[];
for i=1:length(acs)
    if any(~ismember(acs(i)+kernelY, Slines))
        continue;
    end
   acs2(end+1)=acs(i);
end
acs=acs2;
        %ESTIMATE GRAPPA COEFFS---------------------------------
if coef==0
QA = zeros(length(acs)*My, Ncoils*lx*ly); %asumming 3x2 kernel
coef = zeros(lx*ly*Ncoils,Ncoils);
for ii = 1:length(acs),
   ind1 = acs(ii) + kernelY;

   %% determine the k-space filling coefficients 
   for jj=1:length(kernelX),
      y0 = mod( kernelX(jj) + [1:My],My);
      y0( y0 == 0 ) = My;
      QA( My*(ii-1) + (1:My), jj:length(kernelX):size(QA,2) ) = ...
            reshape( permute( Sk(ind1,y0,:) ,[ 2 3 1 ]),My, Ncoils*length(ind1) );
   end
end
    
AA = QA'*QA;

%Better with GNU cgsolv(), Copyright 2001 William Scott Hoge (shoge@ece.neu.edu or shoge@ieee.org) 
% taken from p.529 of Golub and Van Loan, "Matrix Computations," 3rd ed.
% S. 10.2 The Conjugate Gradient Method
% S. 10.2.6 Some Practical Details.
%
  %n(:,l) = cgsolv( AA, A'*b, zeros(size(A,2),1), size(A,2) );  
%
%SUBSTITUTION bicg

for l=1:Ncoils
  b(:,l) = vec( squeeze(Sk( acs, :, l )).' );
% [coef(:,l) TT]=bicg(AA,QA'*b);
 coef(:,l)=cgsolve(AA,QA'*b(:,l),0.01,20,false);
 
%TT Dummy variable to avoid warnings
end

if plot_fit
 figure;plot(abs(vec(b)),abs(vec(QA*coef)),'.');
 
coef2=reshape(coef,[lx,Ncoils,ly,Ncoils]);

kern=zeros(size(Sk,1),size(Sk,2),Ncoils,Ncoils);
kern(-kernelY+floor(Mx/2)+1,-kernelX+floor(My/2)+1,:,:)=permute(coef2,[3,1,2,4]);


sz=size(Sk);
clear tmp;
for j=1:Ncoils
for i=1:Ncoils
tmp(:,:,i,j)=conv2(Sk(:,:,i),kern(:,:,i,j));
end

end

tmp=tmp(sz(1)/2+1:sz(1)*3/2,sz(2)/2+1:sz(2)*3/2,:,:);
tmp=sum(tmp,3);
ftmp=ifft2c(squeeze(tmp))./ifft2c(Sk);
figure;
for i=1:size(ftmp,3)
subplot(4,8,i);
imshow(abs(ftmp(:,:,1)),[0.5,1.2]);
end
figure;
for i=1:size(ftmp,3)
subplot(4,8,i);
imshow(angle(ftmp(:,:,1)),[]);
end
end


if ~do_recon
    Rk=[];
    return;
end


end %If COEF==0
%END ESTIMATION --------------------------------------------


%RECONSTRUCTION---------------------------------------------------

for ii=1:length(ImK),  %For each line

%1.- CONTROL OF INDEX
  if find(Slines==ii), continue, end

  %ELSE: if lines no sampled--> Recons 

  ind1 = ii + kernelY;
  
  if true
  ind1(ind1<1)=ind1(ind1<1)+Mx;
  ind1(ind1>Mx)=ind1(ind1>Mx)-Mx;
  else
    if ( ind1(1) < 1) || (ind1(end) > Mx);    
       continue;%do nothing
    end
  end
        tempo = (ImK(ind1)==ones(length(ind1),1));
  
  if sum(tempo)~=length(ind1)
      continue; 
  end;  
   
%2.- RECONSTRUCT K SPACE

  MA = zeros(My,lx*ly*Ncoils);
  for jj=1:lx
     D1 = mod( kernelX(jj) + [1:My] , My );
     D1(D1==0) = My;   
     MA((1:My),jj:lx:lx*ly*Ncoils)=reshape(permute(squeeze(Sk(ind1,D1,:)),[ 2 3 1 ]),My, Ncoils*ly);
  end %jj
   
  for l=1:Ncoils,   
     Sk(ii,:,l)=MA*coef(:,l); %RECONST
  end

end  %ii
%END RECONSTR----------------------------------------------------------------------

Rk=permute(Sk,[2,1,3]);
%Rk = Sk;
%Rx=(k2x(Rk,1));


% cgsolve.m
%
% Solve a symmetric positive definite system Ax = b via conjugate gradients.
%
% Usage: [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose)
%
% A - Either an NxN matrix, or a function handle.
%
% b - N vector
%
% tol - Desired precision.  Algorithm terminates when 
%    norm(Ax-b)/norm(b) < tol .
%
% maxiter - Maximum number of iterations.
%
% verbose - If 0, do not print out progress messages.
%    Default = 1.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose)

if (nargin < 5), verbose = 1; end

implicit = isa(A,'function_handle');

x = zeros(length(b),1);
r = b;
d = r;
delta = r'*r;
delta0 = b'*b;
numiter = 0;
bestx = x;
bestres = sqrt(delta/delta0); 
while ((numiter < maxiter) & (delta > tol^2*delta0))

  % q = A*d
  if (implicit), q = A(d);  else, q = A*d;  end
 
  alpha = delta/(d'*q);
  x = x + alpha*d;
  
  if (mod(numiter+1,50) == 0)
    % r = b - Aux*x
    if (implicit), r = b - A(x);  else, r = b - A*x;  end
  else
    r = r - alpha*q;
  end
  
  deltaold = delta;
  delta = r'*r;
  beta = delta/deltaold;
  d = r + beta*d;
  numiter = numiter + 1;
  if (sqrt(delta/delta0) < bestres)
    bestx = x;
    bestres = sqrt(delta/delta0);
  end    
  
  if ((verbose) & (mod(numiter,50)==0))
    disp(sprintf('cg: Iter = %d, Best residual = %8.3e, Current residual = %8.3e', ...
      numiter, bestres, sqrt(delta/delta0)));
  end
  
end

if (verbose)
  disp(sprintf('cg: Iterations = %d, best residual = %14.8e', numiter, bestres));
end
x = bestx;
res = bestres;
iter = numiter;




