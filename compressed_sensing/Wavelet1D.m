function res = Wavelet1D(x,level)
% res = Wavelet_ISU(N12,level)
%
% implements a wavelet operator
%N12 is the size of the 2D/3D images.
if ~exist('level','var')
    level=2;
end
dwtmode('per');

N=length(x);
 Xtmp=eye(N);
 W=zeros(N);
  for p=1:N
      W(:,p)=wavedec(Xtmp(:,p),level,'db4');
  end
  
  
  res=W*x(:);
  

  

